use std::ops::Range;

use crate::error::{FastqError, FastqPosition, Result};

#[derive(Clone, Copy)]
pub(crate) struct Line<'a> {
    pub(crate) bytes: &'a [u8],
    pub(crate) start: usize,
}

impl Line<'_> {
    #[inline]
    pub(crate) fn len(self) -> usize {
        self.bytes.len()
    }

    #[inline]
    pub(crate) fn range(self) -> Range<usize> {
        self.start..self.start + self.bytes.len()
    }
}

#[derive(Clone, Copy)]
pub(crate) struct RecordLines<'a> {
    pub(crate) name: Line<'a>,
    pub(crate) seq: Line<'a>,
    pub(crate) plus: Line<'a>,
    pub(crate) qual: Line<'a>,
}

#[derive(Clone, Copy)]
pub(crate) struct RecordValidation {
    require_nonempty_id: bool,
}

impl RecordValidation {
    pub(crate) const DEFAULT: Self = Self {
        require_nonempty_id: true,
    };
    pub(crate) const TRUSTED_PACK: Self = Self {
        require_nonempty_id: false,
    };
}

#[derive(Clone, Copy)]
pub(crate) struct SlabLineLayout {
    pub(crate) line_count: usize,
    pub(crate) complete_lines: usize,
    pub(crate) has_partial_trailing_line: bool,
}

pub(crate) fn slab_line_layout(
    bytes: &[u8],
    newline_offsets: &[usize],
    eof: bool,
) -> SlabLineLayout {
    let has_partial_trailing_line = !eof
        && newline_offsets
            .last()
            .map_or(!bytes.is_empty(), |&nl| nl + 1 < bytes.len());
    let has_final_line = eof
        && newline_offsets
            .last()
            .map_or(!bytes.is_empty(), |&nl| nl + 1 < bytes.len());
    let line_count = newline_offsets.len() + usize::from(has_final_line);
    SlabLineLayout {
        line_count,
        complete_lines: (line_count / 4) * 4,
        has_partial_trailing_line,
    }
}

#[inline]
pub(crate) fn line_start(newline_offsets: &[usize], line: usize) -> usize {
    if line == 0 {
        0
    } else {
        newline_offsets[line - 1] + 1
    }
}

#[inline]
pub(crate) fn line_bounds(bytes: &[u8], newline_offsets: &[usize], line: usize) -> (usize, usize) {
    let start = line_start(newline_offsets, line);
    let end = if line < newline_offsets.len() {
        newline_offsets[line]
    } else {
        bytes.len()
    };
    (start, trim_cr_end(bytes, start, end))
}

#[inline]
pub(crate) fn line<'a>(bytes: &'a [u8], newline_offsets: &[usize], line_index: usize) -> Line<'a> {
    let (start, end) = line_bounds(bytes, newline_offsets, line_index);
    Line {
        bytes: &bytes[start..end],
        start,
    }
}

#[inline]
pub(crate) fn record_lines<'a>(
    bytes: &'a [u8],
    newline_offsets: &[usize],
    line_index: usize,
) -> RecordLines<'a> {
    RecordLines {
        name: line(bytes, newline_offsets, line_index),
        seq: line(bytes, newline_offsets, line_index + 1),
        plus: line(bytes, newline_offsets, line_index + 2),
        qual: line(bytes, newline_offsets, line_index + 3),
    }
}

pub(crate) fn direct_line<'a>(bytes: &'a [u8], cursor: &mut usize, eof: bool) -> Option<Line<'a>> {
    let start = *cursor;
    if start >= bytes.len() {
        return None;
    }

    let mut end = start;
    while end < bytes.len() && bytes[end] != b'\n' {
        end += 1;
    }
    if end == bytes.len() && !eof {
        return None;
    }

    *cursor = if end < bytes.len() { end + 1 } else { end };
    let end = trim_cr_end(bytes, start, end);
    Some(Line {
        bytes: &bytes[start..end],
        start,
    })
}

pub(crate) fn validate_record(
    record: RecordLines<'_>,
    base_offset: u64,
    record_index: u64,
    validation: RecordValidation,
) -> Result<()> {
    if record.name.bytes.first() != Some(&b'@') {
        return Err(format_at(
            "header must start with `@`",
            base_offset,
            record.name.start,
            record_index,
            0,
        ));
    }
    if validation.require_nonempty_id
        && (record.name.len() == 1 || fastq_id_token(record.name.bytes).is_empty())
    {
        return Err(format_at(
            "empty FASTQ id",
            base_offset,
            record.name.start,
            record_index,
            0,
        ));
    }
    if record.plus.bytes.first() != Some(&b'+') {
        return Err(format_at(
            "plus line must start with `+`",
            base_offset,
            record.plus.start,
            record_index,
            2,
        ));
    }
    let seq_len = record.seq.len();
    let qual_len = record.qual.len();
    if seq_len != qual_len {
        return Err(format_at(
            format!("quality length {qual_len} != sequence length {seq_len}"),
            base_offset,
            record.qual.start,
            record_index,
            3,
        ));
    }
    Ok(())
}

fn fastq_id_token(name: &[u8]) -> &[u8] {
    let name = name.strip_prefix(b"@").unwrap_or(name);
    let end = name
        .iter()
        .position(u8::is_ascii_whitespace)
        .unwrap_or(name.len());
    &name[..end]
}

pub(crate) fn fast_slash_pair_ids_match(first_name: &[u8], second_name: &[u8]) -> Option<bool> {
    let first = first_name.strip_prefix(b"@").unwrap_or(first_name);
    let second = second_name.strip_prefix(b"@").unwrap_or(second_name);
    if first.len() >= 3 && first.len() == second.len() && first.ends_with(b"/1") {
        return Some(
            second.ends_with(b"/2") && first[..first.len() - 2] == second[..second.len() - 2],
        );
    }

    let first_end = token_end(first);
    let second_end = token_end(second);
    let first = &first[..first_end];
    let second = &second[..second_end];
    if first.len() < 3 || second.len() < 3 {
        return None;
    }
    let first_suffix = &first[first.len() - 2..];
    let second_suffix = &second[second.len() - 2..];
    if first_suffix != b"/1" || first.len() != second.len() {
        return None;
    }
    Some(second_suffix == b"/2" && first[..first.len() - 2] == second[..second.len() - 2])
}

pub(crate) fn normalized_pair_id(name: &[u8]) -> &[u8] {
    let name = name.strip_prefix(b"@").unwrap_or(name);
    let token = &name[..token_end(name)];
    strip_pair_suffix(token)
}

pub(crate) fn strip_pair_suffix(id: &[u8]) -> &[u8] {
    if id.len() >= 2 && (id.ends_with(b"/1") || id.ends_with(b"/2")) {
        &id[..id.len() - 2]
    } else {
        id
    }
}

#[inline]
pub(crate) fn format_at(
    message: impl Into<String>,
    base_offset: u64,
    local_offset: usize,
    record_index: u64,
    line_index: u8,
) -> FastqError {
    FastqError::FormatAt {
        message: message.into(),
        position: FastqPosition::new(base_offset + local_offset as u64, record_index, line_index),
    }
}

#[inline]
pub(crate) fn trim_cr_end(bytes: &[u8], start: usize, end: usize) -> usize {
    if end > start && bytes[end - 1] == b'\r' {
        end - 1
    } else {
        end
    }
}

#[inline]
fn token_end(bytes: &[u8]) -> usize {
    let mut end = 0;
    while end < bytes.len() && !bytes[end].is_ascii_whitespace() {
        end += 1;
    }
    end
}
