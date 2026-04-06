def get_io_reads(bam_reads, pod5_fh, reverse_signal=False, missing_ok=False):
      pod5_reads = get_pod5_reads(
          pod5_fh, [bam_read.query_name for bam_read in bam_reads]
      )
      io_reads = []
      skipped_count = 0
      for bam_read in bam_reads:
          try:
              if bam_read.query_name not in pod5_reads:
                  if missing_ok:
                      skipped_count += 1
                      continue
                  else:
                      raise RemoraError(
                          f"BAM record '{bam_read.query_name}' not found in
  POD5"
                      )
              io_read = Read.from_pod5_and_alignment(
                  pod5_read_record=pod5_reads[bam_read.query_name],
                  alignment_record=bam_read,
                  reverse_signal=reverse_signal,
              )
          except KeyError:
              if missing_ok:
                  skipped_count += 1
                  continue
              else:
                  raise RemoraError(
                      f"BAM record '{bam_read.query_name}' not found in POD5"
                  )
          except RemoraError as e:
              error_msg = str(e)
              if "Move table" in error_msg and "discordant" in error_msg:
                  skipped_count += 1
                  continue
              if missing_ok:
                  skipped_count += 1
                  continue
              else:
                  raise
          except Exception as e:
              if missing_ok:
                  skipped_count += 1
                  continue
              else:
                  raise
          io_reads.append(io_read)
      return io_reads, skipped_count

