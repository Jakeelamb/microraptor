use microraptor::Result;
use microraptor::benchutil::StreamStats;

use super::{Config, Measurement, measure};

pub(super) fn measurements(input: &[u8], config: &Config) -> Result<Vec<Measurement>> {
    Ok(vec![
        measure_fasta_index_build("build_fasta_index", input, config)?,
        measure_fasta_fetch_repeated("fetch_repeated_range", input, config)?,
        measure_fasta_partition_plan("plan_fasta_partitions", input, config)?,
        measure_fasta_reference_chunks("reference_chunks", input, config)?,
    ])
}

fn measure_fasta_index_build(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    measure(name, input.len(), config.iters, || {
        let index = microraptor::build_fasta_index(std::io::Cursor::new(input))?;
        Ok(index_stats(&index))
    })
}

fn measure_fasta_fetch_repeated(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    let index = microraptor::build_fasta_index(std::io::Cursor::new(input))?;
    let Some(entry) = index.entries().first() else {
        return measure(name, input.len(), config.iters, || {
            Ok(StreamStats::default())
        });
    };
    let range_end = entry.len.min(config.read_len.max(1) as u64);
    let reference_name = entry.name.clone();
    let repeats = config.records.max(1).min(1_024);
    measure(name, input.len(), config.iters, || {
        let mut reader =
            microraptor::IndexedFastaReader::new(std::io::Cursor::new(input), index.clone());
        let mut stats = StreamStats::default();
        for _ in 0..repeats {
            let seq = reader.fetch(&reference_name, 0..range_end)?;
            stats.observe_sequence_record(&reference_name, &seq);
        }
        Ok(stats)
    })
}

fn measure_fasta_partition_plan(name: &str, input: &[u8], config: &Config) -> Result<Measurement> {
    let index = microraptor::build_fasta_index(std::io::Cursor::new(input))?;
    measure(name, input.len(), config.iters, || {
        let partitions = microraptor::plan_fasta_partitions(
            &index,
            microraptor::FastaPartitionConfig::new(config.workers, 31),
        )?;
        let mut stats = StreamStats::default();
        for partition in partitions {
            stats.observe_sequence_record(&partition.name, &[b'N']);
            stats.bases = stats
                .bases
                .checked_add(partition.core_len().saturating_sub(1))
                .ok_or_else(|| {
                    microraptor::FastqError::Format(
                        "partition benchmark base count overflow".into(),
                    )
                })?;
        }
        Ok(stats)
    })
}

fn measure_fasta_reference_chunks(
    name: &str,
    input: &[u8],
    config: &Config,
) -> Result<Measurement> {
    let index = microraptor::build_fasta_index(std::io::Cursor::new(input))?;
    let Some(entry) = index.entries().first() else {
        return measure(name, input.len(), config.iters, || {
            Ok(StreamStats::default())
        });
    };
    let reference_name = entry.name.clone();
    let range = 0..entry.len;
    let chunk_bases = config.read_len.max(1) as u64;
    measure(name, input.len(), config.iters, || {
        let mut reader =
            microraptor::IndexedFastaReader::new(std::io::Cursor::new(input), index.clone());
        let mut stats = StreamStats::default();
        for chunk in reader.reference_chunks(&reference_name, range.clone(), chunk_bases)? {
            let chunk = chunk?;
            stats.observe_sequence_record(&chunk.name, &chunk.seq);
        }
        Ok(stats)
    })
}

fn index_stats(index: &microraptor::FastaIndex) -> StreamStats {
    let mut stats = StreamStats::default();
    for entry in index.entries() {
        stats.observe_sequence_record(&entry.name, &[b'N']);
        stats.bases += entry.len.saturating_sub(1);
    }
    stats
}
