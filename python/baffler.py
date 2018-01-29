"""
Calculate the breakpoint allele frequency (BAF) of structural rearrangements
detected by Delly.
"""

import os
import datetime
import pandas as pd
import numpy as np
import pysam
import SmithWater
import svHelpers

class baffler():

    """
    For each input structural rearrangement, query the BAM file for supporting
    reads, verify the listed breakpoint, and calculate BAF and breakpoint coverage.
    """

    def __init__(self, reference_fasta, log=None, log_level=None):

        """
        Baffler parameters can be configured using the svHelpers.configure
        method. Defaults are set here.
        """

        # ========== Baffler Parameters ========================================

        self.PYSAM_ZERO_BASED = 1

        self.sequence_alignment = {
            'gapopen': 10.0,
            'gapextend': 1.0,
            'identity_percent': 80.0,
            'similarity_percent': 80.0,
            'length_cutoff': 5
        }

        self.ratio_filter = {
            'clipped_reads': 0.8,
            'aligned_reads': 0.8
        }
        self.BAF = {
            'extend_soft_clipped_side': 20,
            'extend_mapped_side': 50,
            'max_junction_diff': 10,
        }
        self.max_coverage = {
            'coverage': 500,
            'window': 1000
        }

        # ======================================================================

        self.baffler_columns = ['passed_ratio_filter',
                                'reads_supporting_translocation',
                                'aligned_breakpoint_1',
                                'clip_pos_breakpoint_1',
                                'COV_TYPE_1', 'BAF_1',
                                'AVG_COV_1', 'COV_STR_1',
                                'aligned_breakpoint_2',
                                'clip_pos_breakpoint_2',
                                'COV_TYPE_2', 'BAF_2',
                                'AVG_COV_2', 'COV_STR_2']

        self.reference_fasta = reference_fasta

        # CIGAR string: see https://samtools.github.io/hts-specs/SAMv1.pdf
        # page 5 for more information.
        pysam_cigar_enums = [[0, 'M'], [1, 'I'], [2, 'D'], [3, 'N'], [4, 'S'],
                             [5, 'H'], [6, 'P'], [7, '='], [8, 'X']]

        self.pysam_cigar_enums = pd.DataFrame.from_records(
            pysam_cigar_enums,
            columns=['operation_enum', 'operation_str']
        )

        # Set up log file
        if log is None:
            log = ('./baffler.tmp.' +
                   str(datetime.datetime.now()).replace(' ', '') +
                   '.log')
        if log_level is None:
            log_level = 'DEBUG'

        log_id = 'baffler.' + str(np.random.randint(0, 1000))
        self.logger = svHelpers.setup_log(log_id=log_id,
                                          logfile=log,
                                          level=log_level)

    def run_baffler(self, tumor_tab, tumor_bam, merge=False, output_tab=None):

        """
        Iterate through input Delly rearrangements and run baffler.baffle to
        calculate BAF.

        Usage:

            baffler.run_baffler(tumor_tab, tumor_bam, merge=False,
                                output_tab=None)

        Arguments:

            tumor_tab (str or pd.DataFrame): either path to Delly rearrangements
                tab file or DataFrame of Delly rearrangements.

            tumor_bam (str): path to BAM file for read lookup.

            merge (bool): Output index (chromosome, position, CT) and
                baffler.baffler_columns with original tumor_tab columns.

            output_tab (str): path to write results tab file.

        Returns:

            results (pd.DataFrame): tumor_tab DataFrame with BAF, AVG_COV,
                and breakpoint classification added, as well as aligned and
                clipped breakpoints for each Delly breakpoint.

        """
        # Create pysam.AlignmentFile object for split read lookup
        tb_obj = pysam.AlignmentFile(tumor_bam, 'rb')
        # Create smithwater object for sequence comparison
        smithwater_obj = SmithWater.SmithWater(self.reference_fasta)

        # ========== Argument Checking =========================================

        if not isinstance(tumor_tab, pd.DataFrame) and not isinstance(tumor_tab, str):
            error_msg = "tumor_tab argument must either be a string or a \
                         Pandas DataFrame"
            raise ValueError(error_msg)

        # ========== Data Extraction ===========================================

        if isinstance(tumor_tab, pd.DataFrame):
            df_tumor = tumor_tab.copy()
        else:
            if not os.path.exists(tumor_tab):
                error_msg = "tumor_tab must be a valid path or a Pandas \
                             DataFrame."
                raise ValueError(error_msg)
            else:
                df_tumor = pd.DataFrame.from_csv(tumor_tab, sep="\t", index_col=False)

        # ========== Data Cleansing ============================================

        df_tumor.Chromosome1 = svHelpers.clean_chromosomes(df_tumor.Chromosome1)
        df_tumor.Chromosome2 = svHelpers.clean_chromosomes(df_tumor.Chromosome2)

        # ========== Calculate BAF =============================================

        df_tumor = self.flag_excessive_coverage(df_tumor, tb_obj)
        temp = []
        for idx, x in df_tumor.iterrows():
            placeholder = x.to_frame().transpose().copy()
            if x.excessive_coverage == True:
                # Excessive coverage - Log as warning to be debugged
                ec_msg = "Excessive coverage found for record index: {idx}, \
                          {chr1}:{pos1};{chr2}:{pos2}"
                ec_msg = ec_msg.format(idx=idx,
                                       chr1=x.Chromosome1,
                                       pos1=x.Position1,
                                       chr2=x.Chromosome2,
                                       pos2=x.Position2)
                self.logger.debug(ec_msg)
                placeholder = self._fill_placeholder(placeholder, 'EC')
                temp.append(placeholder)
                continue
            # Try to calculate the translocation baf, otherwise fill the
            # dataframe with a placeholder
            try:
                result = self.baffle(x.Chromosome1,
                                     int(x.Position1),
                                     x.Chromosome2,
                                     int(x.Position2),
                                     x.CT,
                                     tb_obj,
                                     smithwater_obj)
                if result is not None:
                    temp.append(result)
                else:
                    # Result will be None when:
                    #   - soft-clips could not be aligned via SmithWater
                    #   - resolved breakpoints != 2
                    placeholder = self._fill_placeholder(placeholder, np.nan)
                    temp.append(placeholder)
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                # All unexpected errors should be logged and trigger
                # placeholders.
                error_msg = "Could not calculate translocation BAF for \
                             record index: {idx}, {chr1}:{pos1};{chr2}:{pos2}"
                error_msg = error_msg.format(idx=idx,
                                             chr1=x.Chromosome1,
                                             pos1=x.Position1,
                                             chr2=x.Chromosome2,
                                             pos2=x.Position2)
                self.logger.exception(error_msg)
                placeholder = self._fill_placeholder(placeholder, 'ERROR')
                temp.append(placeholder)

        # ========== Prepare Results ===========================================

        results = pd.concat(temp)
        results.Chromosome1 = svHelpers.clean_chromosomes(results.Chromosome1)
        results.Chromosome2 = svHelpers.clean_chromosomes(results.Chromosome2)

        index_cols = ['Chromosome1', 'Position1', 'Chromosome2', 'Position2', 'CT']
        results = results[index_cols + self.baffler_columns]

        if merge is True:
            results = pd.merge(df_tumor, results, how='left', on=index_cols)

        if output_tab is not None:
            results.to_csv(output_tab, sep="\t", index=False)

        return results

    def _fill_placeholder(self, placeholder, comment):

        """
        Fill baffler columns in placeholder records with comments describing
        reason for error preventing successful BAF calculation.

        Possible errors:

            - EC: At least one of the Delly breakpoints was in a region of
                  excessive coverage

            - ERROR: baffler.baffle raised an unexpected error.

            - None (np.nan): baffler.baffle returned None. Possible reasons:
                * soft clips not aligned via EMBOSS Water
                * more/less than two breakpoints resolved
                See log file for more details.

        Arguments:

            placeholder: placeholder record to fill

            comment: comment describing error, to fill baffler columns with

        Returns:

            placeholder: placeholder record with comment inserted into all
                baffler columns.

        """

        for column in self.baffler_columns:
            if ((column == 'BAF_1') or (column == 'BAF_2')):
                # Set BAF columns to NA so that they can be coerced to float
                placeholder[column] = np.nan
            else:
                placeholder[column] = comment
        return placeholder

    def flag_excessive_coverage(self, df_tumor, tb_obj):

        """
        Flag records where either junction has coverage greater than
        self.max_coverage['coverage'] within (self.max_coverage['window']/2)
        of either side of the breakpoint.

        Usage:

            baffler.flag_excessive_coverage(df_tumor, tb_obj)

        Arguments:

            df_tumor (pd.DataFrame): Rearrangements DataFrame.

            tb_obj (pysam.AlignmentFile): Tumor BAM file.

        Returns:

            df_tumor (pd.DataFrame): Rearrangements DataFrame with column
                'excessive_coverage' added.

        """

        offset = self.max_coverage['window'] / 2
        max_reads = self.max_coverage['coverage']

        df_tumor['junction1_coverage'] = df_tumor.apply(
            lambda x: self._fetch_coverage(x.Chromosome1,
                                           int(x.Position1 - offset),
                                           int(x.Position1 + offset),
                                           tb_obj)[1],
            axis=1)

        df_tumor['junction2_coverage'] = df_tumor.apply(
            lambda x: self._fetch_coverage(x.Chromosome2,
                                           int(x.Position2 - offset),
                                           int(x.Position2 + offset),
                                           tb_obj)[1],
            axis=1)

        df_tumor['excessive_coverage'] = False
        df_tumor.ix[((df_tumor.junction1_coverage > max_reads) |
                     (df_tumor.junction2_coverage > max_reads)),
                    'excessive_coverage'] = True

        return df_tumor

    def baffle(self, chromosome1, position1, chromosome2, position2, CT,
               bam_obj, smithwater_obj):

        """
        Fetch soft-clipped reads on either side of the rearrangement, and align
        their soft-clipped sections to the other side. Based on the clip and
        aligned breakpoints, classify each of the breakpoints as either
        TIER1 or TIER2 and calculate their BAFs.

            clip breakpoint: position that the read is soft-clipped
            aligned breakpoint: position that soft-clip aligns to the other side
                                of the rearrangement

        Returns None when:
            - soft-clips could not be aligned via SmithWater
            - More or less than 2 breakpoints are resolved for the rearrangement

        Usage:

            baffler.baffle(chromosome1, position1, chromosome2, position2, CT,
                           bam_obj, smithwater_obj)

        Arguments:

            chromosome1 (str): Chromosome of first Delly breakpoint.

            position1 (int): Position of first Delly breakpoint.

            chromosome2 (str): Chromosome of second Delly breakpoint.

            position2 (int): Position of second Delly breakpoint.

            CT (str): Delly CT string indicating event connection type. Contains
                the orientation of DNA on each side of rearrangement, delimited
                by 'to'. Example: '3to5'.

            bam_obj (pysam.AlignmentFile): BAM file for read lookup.

            smithwater_obj (SmithWater.SmithWater): SmithWater object for
                sequence alignment via EMBOSS Water.

        Returns:

            results (pd.DataFrame): DataFrame with the following columns:
                'Chromosome1', 'Position1', 'Chromosome2', 'Position2', 'CT' and
                all columns in self.baffler_columns.
                See baffler.merge_junction_output for more information.

        """

        # position1 = int(position1)
        # position2 = int(position2)
        assert isinstance(position1, int), "position1 is not an integer: %r" % position1
        assert isinstance(position2, int), "position2 is not an integer: %r" % position2
        assert isinstance(chromosome1, str), "chromosome1 is not a string: %r" % chromosome1
        assert isinstance(chromosome2, str), "chromosome2 is not a string: %r" % chromosome2
        assert isinstance(CT, str), "CT is not a string: %r" % CT
        assert 'to' in CT, '\'to\' is not in CT: %r' % CT

        or1, or2 = CT.split('to')

        reads_first = self.map_soft_clipped_to_reference(
            CT=CT,
            orientation=or1,
            chromosome=chromosome1,
            position=position1,
            junction_window=10,
            ref_chromosome=chromosome2,
            ref_position=position2,
            ref_window=self.build_ref_window(or2),
            bam_obj=bam_obj,
            smithwater_obj=smithwater_obj
        )

        reads_second = self.map_soft_clipped_to_reference(
            CT=CT,
            orientation=or2,
            chromosome=chromosome2,
            position=position2,
            junction_window=10,
            ref_chromosome=chromosome1,
            ref_position=position1,
            ref_window=self.build_ref_window(or1),
            bam_obj=bam_obj,
            smithwater_obj=smithwater_obj
        )

        if (len(reads_first) == 0 and len(reads_second) == 0):
            error_msg = "No split reads found mapping to other side \
                         on either junction for: \
                        {chr1}:{pos1};{chr2}:{pos2}. Skipping..."
            error_msg = error_msg.format(chr1=chromosome1, pos1=position1,
                                         chr2=chromosome2, pos2=position2)
            self.logger.debug(error_msg)

        # Check that at least one side of the rearrangement has soft clips that
        # were aligned via SmithWater.
        if ('identity_percent' not in reads_first.columns and
                'identity_percent' not in reads_second.columns):
            error_msg = "SmithWater Alignment columns not found for reads on \
                         either junction. Error in alignment to reference for: \
                        {chr1}:{pos1};{chr2}:{pos2}. Skipping..."
            error_msg = error_msg.format(chr1=chromosome1, pos1=position1,
                                         chr2=chromosome2, pos2=position2)
            self.logger.debug(error_msg)
            return None

        if len(reads_first.columns) > len(reads_second.columns):
            # Second junction has no soft clipped reads, but make sure
            # reads_second has columns as if it did
            reads_second = reads_second.reindex_axis(reads_first.columns, axis=1)

        elif len(reads_first.columns) < len(reads_second.columns):
            # First junction has no soft clipped reads, but make sure
            # reads_first has columns as if it did
            reads_first = reads_first.reindex_axis(reads_second.columns, axis=1)

        # Add a column to designate which side of the translocation Delly calls
        # '1' and which side Delly calls '2' so we can decouple them after we
        # put them together.

        reads_first['event_order'] = 1
        reads_second['event_order'] = 2

        sc_first = reads_first.pipe(self.apply_mapping_filters)
        sc_second = reads_second.pipe(self.apply_mapping_filters)
        
        if len(sc_first) < 1 and len(sc_second) < 1:
            error_msg = "All soft-clipped reads removed by mapping filters for: \
                        {chr1}:{pos1};{chr2}:{pos2}. Skipping..."
            error_msg = error_msg.format(chr1=chromosome1, pos1=position1,
                                         chr2=chromosome2, pos2=position2)
            self.logger.debug(error_msg)
            return None
        elif len(sc_first) < 1:
            sc_first = reads_first.reindex_axis(sc_second.columns, axis=1)
        elif len(sc_second) < 1:
            sc_second = reads_second.reindex_axis(sc_first.columns, axis=1)

        sc_first = self.mark_aligned_breakpoints(sc_first, opp_orientation=or2)
        sc_second = self.mark_aligned_breakpoints(sc_second, opp_orientation=or1)

        sc_reads = pd.concat([sc_first, sc_second])
        breakpoints = self.resolve_breakpoints(sc_reads)
        breakpoints = breakpoints.sort(columns=['event_order'],
                                       ascending=True)

        if len(breakpoints) != 2:
            # While this should never happen, the next section depends on there
            # being only two breakpoints (one on either side of the
            # rearrangement), so we check explicitly.
            error_msg = "More or less than two junctions resolved for \
                         {chr1}:{pos1};{chr2}:{pos2}. Skipping..."
            error_msg = error_msg.format(chr1=chromosome1, pos1=position1,
                                         chr2=chromosome2, pos2=position2)
            self.logger.debug(error_msg)
            return None

        results = (breakpoints
                   .pipe(self.swap_aligned_breakpoints)
                   .pipe(self.calculate_baf, sc_reads, bam_obj)
                   .pipe(self.merge_junction_output)
                   .pipe(self.apply_ratio_filter))

        results['CT'] = CT

        return results

    def map_soft_clipped_to_reference(self, CT, orientation, chromosome,
                                      position, junction_window, ref_chromosome,
                                      ref_position, ref_window, bam_obj,
                                      smithwater_obj):
        """
        Fetch all soft-clipped reads within 'junction_window' bases of the input
        'chromosome' and 'position' and align their soft-clips to the input
        region of the reference genome. 'chromosome' and 'position' correspond
        to 'this' side of the rearrangement, while 'ref_chromosome' and
        'ref_position' correspond to the 'other' side of the rearrangement.

        See the SmithWater module for more information on pairwise sequence
        alignment using EMBOSS Water.

        Usage:

            baffler.map_soft_clipped_to_reference(CT, orientation, chromosome,
                                                  position, junction_window,
                                                  ref_chromosome, ref_position,
                                                  ref_window, bam_obj,
                                                  smithwater_obj)

        Arguments:

            CT (str): Delly CT string indicating event connection type. Contains
                the orientation of DNA on each side of rearrangement, delimited
                by 'to'. Example: '3to5'.

            chromosome (str): chromosome to search for soft-clipped reads

            position (int): position to search for soft-clipped reads

            junction_window (int): distance to search around position for
                soft-clipped reads

            ref_chromosome (str): chromosome to align soft-clips

            ref_position (int): position to align soft-clips

            ref_window ((int, int)): tuple of distances (left, right) around
                ref_position to align soft-clips. Passed as tuple to allow uneven
                window around ref_position, depending on junction orientation.

            orientation (str): Orientation of the DNA on 'this' side of the
                rearrangement. Either '5' (forward) or '3' (reverse).

            bam_obj (pysam.AlignmentFile): BAM file for read lookup.

            smithwater_obj (SmithWater.SmithWater): SmithWater object for
                sequence alignment via EMBOSS Water.


        Returns:

            soft_clipped (pd.DataFrame): soft-clipped reads and their alignment
                to the specified region of the reference genome.

        """

        fetch_start = position - junction_window
        fetch_end = position + junction_window

        ref_window_left, ref_window_right = ref_window
        reference_start = ref_position - ref_window_left
        reference_end = ref_position + ref_window_right

        soft_clipped = self.get_soft_clipped_reads(orientation=orientation,
                                                   chromosome=chromosome,
                                                   position=position,
                                                   start=fetch_start,
                                                   stop=fetch_end,
                                                   bam_obj=bam_obj)

        prime1, prime2 = CT.split('to')

        if prime1 == prime2:
            revcomp = (True, True)
        else:
            revcomp = (False, False)

        alignments = []
        for idx, row in soft_clipped.iterrows():
            if pd.isnull(row.seq): continue
            alignments.append(
                smithwater_obj.align_to_reference(
                    query_name=row.query_name,
                    sequence=row.seq,
                    reference_chromosome=ref_chromosome,
                    reference_start=reference_start,
                    reference_end=reference_end,
                    reverse=revcomp[0],
                    complement=revcomp[1],
                    gapopen=self.sequence_alignment['gapopen'],
                    gapextend=self.sequence_alignment['gapextend']
                )
            )

        if len(alignments) > 0:
            alignments_df = pd.concat(alignments)
            soft_clipped = pd.merge(soft_clipped,
                                    alignments_df,
                                    how='left',
                                    on='query_name')

        return soft_clipped

    def build_ref_window(self, orientation):

        """
        When aligning soft-clips to the reference genome on the 'other' side
        of the rearrangement, we want to extend self.BAF['extend_soft_clipped_side']
        bases into the soft-clipped region and self.BAF['extend_mapped_side']
        into the mapped region around the breakpoint.

        Since the order of the mapped / soft-clipped side depends on
        orientation, we need a dynamically sized window.

        Usage:

            baffler.build_ref_window(orientation)

        Arguments:

            orientation (str): Orientation of the DNA on the 'other' side of the
                rearrangement, where we want to align the soft-clipped segments.
                Either '5' (forward) or '3' (reverse).

        Returns:

            ref_window ((int, int)): Tuple of distances to extend (left, right)
                around the 'other' Delly breakpoint when aligning soft-clips to
                the reference.

        """

        if orientation == '5':
            return (self.BAF['extend_soft_clipped_side'],
                    self.BAF['extend_mapped_side'])

        if orientation == '3':
            return (self.BAF['extend_mapped_side'],
                    self.BAF['extend_soft_clipped_side'])

    def get_soft_clipped_reads(self, orientation, chromosome, position, start, stop, bam_obj):

        """
        Fetch soft-clipped reads between 'start' and 'stop' on the specified
        'chromosome' from the 'bam_obj'.

        Only soft-clipped reads that correspond to the Delly event orientation
        are selected:

            - If orientation is '5': CIGAR is of the form '*S*M--'.

            - If orientation is '3': CIGAR is of the form '--*M*S'.

        We compute the position that the read is soft-clipped ourselves
        (rather than using pysam) to ensure that it agrees with the Delly
        orientation. This position, 'clip_pos', is the clip breakpoint.

        Usage:

            baffler.get_soft_clipped_reads(orientation, chromosome, position,
                                           start, stop, bam_obj)

        Arguments:

            chromosome (str): chromosome to search for soft-clipped reads.

            position (int): position to search for soft-clipped reads.

            start (int): start of window to search for soft-clipped reads.

            stop (int): end of window to search for soft-clipped reads.

            bam_obj (pysam.AlignmentFile): BAM file for read lookup.

            orientation (str): Orientation of the DNA on 'this' side of the
                rearrangement. Either '5' (forward) or '3' (reverse).

        Returns:

            soft_clipped (pd.DataFrame): DataFrame of soft-clipped read
                information and the clip breakpoint. Columns: 'chromosome',
                'position', 'orientation', 'query_name', 'cigar', 'seq',
                and 'clip_pos'.

        """

        # if chromosome != 'X' and chromosome != 'Y':
        #     chromosome = int(chromosome)
        assert isinstance(chromosome, str), "chromosome is not a string: %r" % chromosome
        assert isinstance(start, int), "start is not an int: %r" % start
        assert isinstance(stop, int), "stop is not an int: %r" % stop

        sc_temp = []
        for read in (bam_obj.fetch(chromosome,
                                   int(start - self.PYSAM_ZERO_BASED),
                                   int(stop - self.PYSAM_ZERO_BASED))):

            if read.cigartuples is None:
                continue
            if 'S' not in read.cigarstring:
                continue

            # pysam.calignedsegment.AlignedSegment.cigartuples returns a list of
            # tuples. Each tuple has form (operation, length). The operation is
            # given as an enum ranging from 0-8. We call this enum the
            # operation_enum. Each enum is associated with a standard CIGAR
            # letter, that we call operation_str, such as 0:'M' for match. The
            # length is the number of bases with that alignment operation.
            # See http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
            # for more information.

            cigar = pd.DataFrame(read.cigartuples,
                                 columns=['operation_enum', 'length'])

            # Get the operation_str for each operation_enum. For example:
            # operation_enum: 0 --> operation_str: 'M' --> Match
            cigar = pd.merge(cigar,
                             self.pysam_cigar_enums,
                             how='left',
                             on='operation_enum')

            clipped_cigar = (cigar
                             .pipe(self.vectorize_cigar, read)            # Vectorize CIGAR (to length of read)
                             .pipe(self.vectorize_sequence, read)         # Vectorize corresponding sequence
                             .pipe(self.clip_cigar, cigar, orientation))  # Clip vectorized CIGAR at start of soft-clip

            if clipped_cigar is None:
                continue # no soft clip

            clipped_seq = ''.join(clipped_cigar.base)                       # Sequence that is soft-clipped
            clip_pos = self.get_clip_breakpoint(clipped_cigar, orientation) # Clip breakpoint

            sc_temp.append([chromosome, position, orientation, read.query_name,
                            read.cigarstring, clipped_seq, clip_pos])

        if len(sc_temp) < 1:
            sc_temp = [[chromosome, position, orientation, np.nan, np.nan,
                        np.nan, np.nan]]

        soft_clipped = pd.DataFrame(sc_temp, columns=['chromosome', 'position',
                                                      'orientation', 'query_name',
                                                      'cigar', 'seq', 'clip_pos'])
        return soft_clipped

    def vectorize_cigar(self, cigar, read):

        """
        Vectorize the CIGAR string so that each vector element represents a
        CIGAR operation on a particular position of the reference genome.

        Add the reference genome positions corresponding to each CIGAR operation
        as well.

        Note that the CIGAR string is passed as a DataFrame where each row
        represents a different CIGAR operation (similar to pysam's cigartuples).

        The CIGAR string '3M2S4M' would be passed as DataFrame:

            operation_enum | operation_str | length
           -----------------------------------------
                   0       |      M        |    3
                   4       |      S        |    2
                   0       |      M        |    4

        And would be converted to:

                cig: ['M', 'M', 'M', 'S', 'S', 'M', 'M', 'M', 'M']

        If the start position were 154 then the associated positions would be:

                pos: [154, 155, 156, 157, 158, 159, 160, 161, 162]

        The output DataFrame has columns 'cig' and 'pos'.

        Usage:

            baffler.vectorize_cigar(cigar, read)

        Arguments:

            cigar (pd.DataFrame): CIGAR string DataFrame. Similar output format
                as pysam.calignedsegment.AlignedSegment.cigartuples but with
                CIGAR operations as strings as well as enums.

            read (pysam.calignedsegment.AlignedSegment): Read object returned
                by pysam.calignmentfile.AlignmentFile.fetch.

        Returns:

            exploded_cigar (pd.DataFrame): Vectorized CIGAR string with a
                'cig' column indicating the CIGAR operation, on the reference
                genome position indicated by the 'pos' column.

        """

        start = read.reference_start - read.query_alignment_start + self.PYSAM_ZERO_BASED

        positions = []
        operations = []
        for idx, row in cigar.iterrows():
            # CIGAR: [operation_str: M, length: 3] --> [M, M, M]
            operation = [row.operation_str] * row.length
            operations.extend(operation)
            # positions: start = 154 --> [154, 155, 156]
            pos = list(range(start, start + row.length))
            positions.extend(pos)
            start = start + row.length

        exploded_cigar = pd.DataFrame({'pos': positions, 'cig': operations})
        return exploded_cigar

    def vectorize_sequence(self, exploded_cigar, read):

        """
        Add the read sequence to the vectorized CIGAR string.

        If for example, the sequence 'ACTAA' had CIGAR '3M2S', starting from
        position 324, the output DataFrame will be:

            pos | cig | base
           ------------------
            154 |  M  |  A
            155 |  M  |  C
            156 |  M  |  T
            157 |  S  |  A
            158 |  S  |  A

        Hard clip and deletions do not have corresponding sequence and are set
        to NA.

        Usage:

            baffler.vectorize_sequence(exploded_cigar, read)

        Arguments:

            exploded_cigar (pd.DataFrame): Vectorized CIGAR string with a
                'cig' column indicating the CIGAR operation, on the reference
                genome position indicate by the 'pos' column.

            read (pysam.calignedsegment.AlignedSegment): Read object returned
                by pysam.fetch

        Returns:

            exploded_cigar (pd.DataFrame): Vectorized CIGAR string. The column
                'base' indicates the aligned base from the input read for each
                position.
        """

        sequence_iterator = iter(list(read.seq))
        temp = []
        for idx, row in exploded_cigar.iterrows():
            if row.cig == 'H' or row.cig == 'D':
                # Hard clips (H) and deletions (D) won't have any sequence
                temp.append(np.nan)
            else:
                try:
                    base = next(sequence_iterator)
                except StopIteration:
                    # Reached the end of the read
                    break
                temp.append(base)
        exploded_cigar['base'] = temp
        return exploded_cigar

    def clip_cigar(self, exploded_cigar, cigar, orientation):

        """
        Select the section of the vectorized CIGAR string that should map to the
        'other' side of the rearrangement (i.e. the soft clip).

        Usage:

            baffler.clip_cigar(cigar, exploded_cigar, orientation)

        Arguments:

            exploded_cigar (pd.DataFrame): Vectorized CIGAR string with a
                'cig' column indicating the CIGAR operation, on the reference
                genome position indicated by the 'pos' column. The column
                'base' indicates the aligned base from the input read for each
                position.

            cigar (pd.DataFrame): CIGAR string DataFrame. Similar output format
                as pysam.calignedsegment.AlignedSegment.cigartuples but with
                CIGAR operations as strings as well as enums.

            orientation (str): Orientation of the DNA on 'this' side of the
                rearrangement. Either '5' (forward) or '3' (reverse).

        Returns:

            exploded_sc (pd.DataFrame): Section of vectorized CIGAR string that
                represents the soft-clipped portion of the read.

        """

        exploded_cigar = exploded_cigar.copy() # Copy so we do not add cigar_section to original
        cigar = cigar.copy() # Copy so we do not add cigar section to original

        # TODO: explain this logic
        cigar['cigar_section'] = ((cigar.operation_str != cigar.operation_str.shift(-1))).shift(1).fillna(False).cumsum()
        exploded_cigar['cigar_section'] = ((exploded_cigar.cig != exploded_cigar.cig.shift(-1))).shift(1).fillna(False).cumsum()

        if orientation == '5':
            # 5 prime junction should have cigar like 50S100M***
            soft_clip = cigar[((cigar.operation_str == 'S') &
                               (cigar.shift(-1).operation_str == 'M'))]
            if len(soft_clip) > 0:
                # The soft clip portion should be at the beginning of the read
                # in 5-prime splits, so choose the first part of the CIGAR
                # which has 'S' followed by 'M'.
                soft_clip = soft_clip[(soft_clip.cigar_section ==
                                       soft_clip.cigar_section.min())]

        elif orientation == '3':
            # 3 prime junction should have cigar like ***100M50S
            soft_clip = cigar[((cigar.operation_str == 'S') &
                               (cigar.shift(1).operation_str == 'M'))]
            if len(soft_clip) > 0:
                # The soft clip portion should be at the end of the read in
                # 3-prime splits, so choose the last part of the CIGAR
                # which has 'M' followed by 'S'.
                soft_clip = soft_clip[(soft_clip.cigar_section ==
                                       soft_clip.cigar_section.max())]

        if len(soft_clip) != 1:
            # If the read is matched and soft-clipped more than once, we choose
            # to skip it. Similarly, skip if there is no 'match + soft clip'
            # sequence.
            return None
        elif soft_clip.operation_str.item() != 'S':
            # If the soft_clip does not have CIGAR operation 'S' then something
            # broke in the orientation logic and the read should be skipped.
            return None

        soft_clip_section = soft_clip.cigar_section.item()

        # Select the soft-clip from the vectorized CIGAR
        exploded_sc = exploded_cigar[exploded_cigar.cigar_section == soft_clip_section]
        exploded_sc = exploded_sc.drop('cigar_section', axis=1)

        return exploded_sc

    def get_clip_breakpoint(self, exploded_sc, orientation):

        """
        Identify the breakpoint on 'this' side of the rearrangement based on the
        position that the read is soft-clipped.

        Usage:

            baffler.get_clip_breakpoint(exploded_sc, orientation)

        Arguments:

            exploded_sc (pd.DataFrame): Section of vectorized CIGAR string that
                represents the soft-clipped portion of the read.

            orientation (str): Orientation of the DNA on 'this' side of the
                rearrangement. Either '5' (forward) or '3' (reverse).

        Returns:

            clip_pos (int): The rearrangement breakpoint based on the position
                of the soft-clip.

        """

        if orientation == '5':
            # The clip breakpoint is the position after the end of the soft-clip
            # when 'this' side of the rearrangement has 5 prime orientation
            clip_pos = exploded_sc.pos.max() + 1
        elif orientation == '3':
            # The clip breakpoint is the position before the beginning of the
            # soft-clip when 'this' side of the rearrangement has 3 prime
            # orientation
            clip_pos = exploded_sc.pos.min() - 1

        return clip_pos

    def apply_mapping_filters(self, reads):

        """
        Apply filters on soft clip alignment quality. Thresholds are set at the
        class-level and can be specified in the YAML config. Filters with
        thresholds that are 'None' will not be applied. Any read failing a
        filter is removed.

        The three filters are:

            - identity: percentage of identical matches between the soft clip
                        and the reference over the aligned region (including
                        gaps)

            - similarity: percentage of matches between the soft clip and the
                          reference over the aligned region (including gaps)

            - length: the size of the soft clip

        Identity and similarity are measures of alignment quality that are
        returned by EMBOSS Water. For more information, refer to:
        http://embossgui.sourceforge.net/demo/manual/water.html

        Usage:

            baffler.apply_mapping_filters(reads)

        Arguments:

            reads (pd.DataFrame): Soft clipped reads. Clipped segments are
                aligned to the other side of the rearrangement using EMBOSS
                water in SmithWater.

        Returns:

            reads (pd.DataFrame): Filtered soft clipped reads with new columns:
                'identity_filter', 'similarity_filter', 'length_filter'.

        """

        idt_pct = self.sequence_alignment['identity_percent']
        if idt_pct is not None:
            reads['identity_filter'] = False
            reads.ix[(reads.identity_percent > idt_pct),
                     'identity_filter'] = True

        sim_pct = self.sequence_alignment['similarity_percent']
        if sim_pct is not None:
            reads['similarity_filter'] = False
            reads.ix[(reads.identity_percent > sim_pct),
                     'similarity_filter'] = True

        length_cutoff = self.sequence_alignment['length_cutoff']
        if length_cutoff is not None:
            reads['length_filter'] = False
            reads.ix[[len(x) >= length_cutoff for x in reads.sequence if pd.notnull(x)],
                     'length_filter'] = True

        # Select valid soft-clips (i.e. those passing all filters)
        sc = reads.copy()
        sc = sc[((sc.identity_filter == True) &
                 (sc.similarity_filter == True) &
                 (sc.length_filter == True))]

        return sc

    def mark_aligned_breakpoints(self, soft_clipped, opp_orientation):

        """
        The aligned breakpoint is the position on the reference genome that the
        first base in the soft-clip aligns to on the other side of the
        rearrangement.

        The breakpoint is dependent on the orientation of the other side of the
        rearrangement:

            Orientation | Soft Clip Position
           ----------------------------------
            Forward (5) | Start
            Reverse (3) | End

        Usage:

            baffler.mark_aligned_breakpoints(soft_clipped, opp_orientation)

        Arguments:

            soft_clipped (pd.DataFrame): Soft clipped reads. Clipped segments
                are aligned to the other side of the rearrangement using EMBOSS
                water in SmithWater.

            opp_orientation (str): Orientation of the DNA on the opposite side
                of the rearrangement. Either '5' (forward) or '3' (reverse).

        Returns:

            soft_clipped (pd.DataFrame): Soft clipped reads with the column
                'aligned' added to indicate the aligned breakpoint.

        """

        soft_clipped['aligned'] = np.nan
        if opp_orientation == '5':
            soft_clipped['aligned'] = soft_clipped.ref_genomic_start
        elif opp_orientation == '3':
            soft_clipped['aligned'] = soft_clipped.ref_genomic_end
        return soft_clipped


    def resolve_breakpoints(self, soft_clipped):

        """
        Select the mode clipped and aligned breakpoints for either end of the
        rearrangement. Count the number of soft-clipped reads with mode
        breakpoints, as well as the number of reads clipped/aligned differently.

        Also count the number of reads that support the translocation: reads
        that have both mode clipped and aligned breakpoints.

        Usage:

            baffler.resolve_breakpoint(soft_clipped)

        Arguments:

            soft_clipped (pd.DataFrame): Soft clipped reads. Clipped segments
                are aligned to the other side of the rearrangement using EMBOSS
                water in SmithWater.

        Returns:

            mode_breakpoints (pd.DataFrame): The mode breakpoints for both ends
                of the rearrangement, the number of soft-clips with these
                breakpoints, and the number of soft-clips with different
                breakpoints (both clipped and aligned).

        """
        # A (20170420): pd.notnull(x) must be checked in the following lambdas
        # since there are several ways that non-aligned soft clips make it here

        # breakpoint: the mode breakpoint (clip/aligned) in the soft-clips.
        # Note we use value_counts instead of mode since pd.mode() returns None
        # if all x have fewer than 2 observations. We would rather just choose
        # one arbitrarily.
        bmode = lambda x: x.value_counts().index[0] if all(pd.notnull(x)) else None
        # breakpoint_reads: the number of soft clips with the mode breakpoint.
        mode_count = lambda x: x.value_counts().values[0] if all(pd.notnull(x)) else None
        # other: The number of soft clips with non-mode breakpoints.
        non_mode_count = lambda x: x.value_counts().values[1:].sum() if all(pd.notnull(x)) else None

        mode_breakpoints = (soft_clipped
                            .groupby(['chromosome', 'position', 'event_order'])
                            .agg({'clip_pos': {'breakpoint': bmode,
                                               'breakpoint_reads': mode_count,
                                               'other': non_mode_count},

                                  'aligned':  {'breakpoint': bmode,
                                               'breakpoint_reads': mode_count,
                                               'other': non_mode_count},

                                  'query_name': 'count'})
                            .reset_index())

        # Reduce the DataFrame header to one level by joining the two levels
        # with an underscore.
        mode_breakpoints.columns = ['_'.join(x).strip('_') for x in mode_breakpoints.columns]

        # TODO: This is risky because it may happen that the breakpoints on
        # either side of the rearrangement happen to be the same number and this
        # would break.
        sc_both_mode = len(
            soft_clipped[soft_clipped.clip_pos.isin(mode_breakpoints.clip_pos_breakpoint) &
                         soft_clipped.aligned.isin(mode_breakpoints.aligned_breakpoint)]
        )

        mode_breakpoints['reads_supporting_translocation'] = sc_both_mode

        mode_breakpoints = mode_breakpoints.rename(
            columns={
                'query_name_count': 'total_split_reads',
                'clip_pos_other': 'clipped_differently',
                'aligned_other': 'aligned_differently'
            }
        )

        return mode_breakpoints

    def swap_aligned_breakpoints(self, breakpoints):

        """
        Aligned breakpoints indicate where the soft-clipped reads from 'this'
        side of the rearrangement align to the 'other' side. Since we want to
        compare aligned and clip breakpoints on the same chromosome, we swap
        the aligned breakpoints.

        Usage:

            baffler.swap_aligned_breakpoints(breakpoints)

        Arguments:

            breakpoints (pd.DataFrame): a 2-row DataFrame containing the
                two Delly breakpoints and their corresponding clip and aligned
                breakpoints.

        Returns:

            breakpoints (pd.DataFrame): breakpoints DataFrame with the column
                'opp_aligned_breakpoint' added.

        """

        breakpoints['opp_aligned_breakpoint'] = list(breakpoints.aligned_breakpoint)[::-1]
        return breakpoints

    def calculate_baf(self, breakpoints, sc_reads, bam_obj):

        """
        Classify Delly breakpoints into two tiers and then calculate BAF
        and breakpoint coverage.

        Delly breakpoints are classified into two tiers based on their
        clip and aligned breakpoints:

        Tier 1:

            - clip breakpoint is not None and;
            - aligned breakpoint is not None and;
            - clip breakpoint and aligned breakpoint are within
              self.BAF['max_junction_diff'] bases of one another

        Tier 2:

            - Any of the above conditions are False

        The numerator of the BAF is the same in both cases: the number of reads,
        passing the alignment filters, which have both mode clip breakpoint
        and mode aligned breakpoint.

        The denominator of the BAF depends on the breakpoint tier.

        Tier 1:

            - The average coverage between the clip breakpoint and the aligned
              breakpoint

        Tier 2:

            - The coverage at the Delly breakpoint

        When fetching coverage, the query names of reads in the BAF numerator
        are passed as baseline_reads so that they are added to the BAF
        denominator regardless of whether they are found by pysam.pileup.

        Usage:

            baffler.calculate_baf(breakpoints, sc_reads, bam_obj)

        Arguments:

            breakpoints (pd.DataFrame): a 2-row DataFrame containing the
                two Delly breakpoints and their corresponding clip and aligned
                breakpoints.

            sc_reads (pd.DataFrame): the soft-clipped reads from either end of
                the rearrangement that pass the alignment filters.

            bam_obj (pysam.AlignmentFile): BAM file for read lookup.

        Returns:

            breakpoints (pd.DataFrame): Delly breakpoints for the rearrangement,
                plus clip and aligned breakpoints, with additional columns:
                'BAF', 'AVG_COV', 'COV_STR', and 'COV_TYPE'.

        """

        # TODO: again, this is risky if the breakpoints happen to be the same
        # on either side of the rearrangement.

        translocation_query_names = sc_reads[
            sc_reads.clip_pos.isin(breakpoints.clip_pos_breakpoint) &
            sc_reads.aligned.isin(breakpoints.aligned_breakpoint)
        ].query_name.tolist()

        for idx, row in breakpoints.iterrows():
            clipped = row.clip_pos_breakpoint
            aligned = row.opp_aligned_breakpoint
            junction_diff = np.abs(clipped - aligned)

            if (pd.isnull(clipped) or pd.isnull(aligned) or
                    (junction_diff > self.BAF['max_junction_diff'])):
                coverage_type = 'TIER2'
                coverage_str, coverage = self._fetch_coverage(
                    row.chromosome,
                    int(row.position),
                    int(row.position+1),
                    bam_obj,
                    baseline_reads=translocation_query_names)
            else:
                coverage_type = 'TIER1'
                start = int(np.min([clipped, aligned]))
                end = int(np.max([clipped, aligned]))
                coverage_str, coverage = self._fetch_coverage(
                    row.chromosome,
                    start,
                    end,
                    bam_obj,
                    baseline_reads=translocation_query_names)

            breakpoints.ix[idx, 'AVG_COV'] = coverage
            breakpoints.ix[idx, 'COV_STR'] = coverage_str
            breakpoints.ix[idx, 'COV_TYPE'] = coverage_type

        breakpoints['BAF'] = np.round(
            (breakpoints.reads_supporting_translocation * 1.0 /
             breakpoints.AVG_COV).values, 4
        )

        return breakpoints

    def merge_junction_output(self, breakpoints):

        """
        Transform 2-row breakpoints DataFrame, used to calculate BAFs, to a
        single record similar to original Delly output.

        Usage:

            baffler.merge_junction_output(breakpoints)

        Arguments:

            breakpoints (pd.DataFrame): a 2-row DataFrame containing the
                two Delly breakpoints and their corresponding clip and aligned
                breakpoints.

        Returns:

            junctions_out (pd.DataFrame): a single-row DataFrame containing
                the Delly breakpoints and baffler columns, identified as '1'
                and '2' based on their original event order in Delly.

        """
        breakpoints['swap_chromosome'] = list(breakpoints.chromosome)[::-1]
        breakpoints['swap_position'] = list(breakpoints.position)[::-1]
        junctions_out = pd.merge(breakpoints,
                                 breakpoints,
                                 how='left',
                                 left_on=['chromosome', 'position'],
                                 right_on=['swap_chromosome', 'swap_position'],
                                 suffixes=['_1', '_2'])

        junctions_out = junctions_out.head(1)
        junctions_out = junctions_out.drop(['event_order_1', 'event_order_2',
                                            'swap_chromosome_1',
                                            'swap_position_1',
                                            'swap_chromosome_2',
                                            'swap_position_2',
                                            'opp_aligned_breakpoint_1',
                                            'opp_aligned_breakpoint_2',
                                            'reads_supporting_translocation_2'],
                                           axis=1)

        junctions_out = junctions_out.rename(
            columns={'chromosome_1': 'Chromosome1',
                     'chromosome_2': 'Chromosome2',
                     'position_1': 'Position1',
                     'position_2': 'Position2',
                     'reads_supporting_translocation_1': 'reads_supporting_translocation'})

        return junctions_out

    def apply_ratio_filter(self, breakpoints):

        """
        Identify whether rearrangement passes ratio filter.

        To pass the ratio filter, both breakpoints must:

            1) Have at least self.ratio_filter['aligned_reads'] percent of
               soft-clips align to the mode aligned breakpoint

            2) Have at least self.ratio_filter['clipped_reads'] percent of
               soft-clips clipped the mode clip breakpoint

        Usage:

            baffler.apply_ratio_filter(breakpoints)

        Arguments:

            breakpoints: Delly breakpoints for the rearrangement,
                plus clip and aligned breakpoints.

        Returns:

            breakpoints: Delly breakpoints for the rearrangement,
                plus clip and aligned breakpoints, with additional columns:
                'passed_ratio_filter'.

        """
        breakpoints['passed_ratio_filter'] = False
        breakpoints.ix[
            ((breakpoints.aligned_breakpoint_reads_1 * 1.0 /
              (breakpoints.aligned_differently_1 + breakpoints.aligned_breakpoint_reads_1)) >=
             self.ratio_filter['aligned_reads']) &

            ((breakpoints.aligned_breakpoint_reads_2 * 1.0 /
              (breakpoints.aligned_differently_2 + breakpoints.aligned_breakpoint_reads_2)) >=
             self.ratio_filter['aligned_reads']) &

            ((breakpoints.clip_pos_breakpoint_reads_1 * 1.0 /
              (breakpoints.clipped_differently_1 + breakpoints.clip_pos_breakpoint_reads_1)) >=
             self.ratio_filter['clipped_reads']) &

            ((breakpoints.clip_pos_breakpoint_reads_2 * 1.0 /
              (breakpoints.clipped_differently_2 + breakpoints.clip_pos_breakpoint_reads_2)) >=
             self.ratio_filter['clipped_reads'])

            , 'passed_ratio_filter'] = True

        return breakpoints

    def _fetch_coverage(self, chromosome, position1, position2, bam_obj,
                        baseline_reads=None):

        """
        Fetch coverage at each base, and average coverage, from the 'bam_obj'
        between 'position1' and 'position2' on 'chromosome' using pysam's pileup
        functionality.

        'baseline_reads' ensure that reads in the BAF numerator always count
        towards the BAF denominator (coverage). The length of the
        list is the baseline coverage at each base between 'position1' and
        'position2'. Any reads found to be covering a base between the two
        positions via pysam pileup are added to the baseline coverage, so
        long as their query_name is not already in 'baseline_reads'.

        Usage:

            baffler._fetch_coverage(chromosome, position1, position2, bam_obj,
                                    baseline_reads=None)

        Arguments:

            chromosome (str): chromosome to fetch coverage

            position1 (int): start position to fetch coverage

            position2 (int): end position to fetch coverage

            bam_obj (pysam.AlignmentFile): BAM file for read lookup.

            baseline_reads (list): List of reads couting towards BAF numerator
                to establish coverage baseline.

        Returns:

            tuple (cov_str, avg_cov):

                cov_str: A string showing coverage at each base in between
                         'position1' and 'position2'

                avg_cov: Average coverage between 'position1' and 'position2'

        """

        # position1 = int(position1)
        # position2 = int(position2)
        assert isinstance(position1, int), "position1 is not an integer: %r" % position1
        assert isinstance(position2, int), "position2 is not an integer: %r" % position2
        assert isinstance(chromosome, str), "chromosome is not a string: %r" % chromosome

        if baseline_reads is None:
            baseline_reads = []

        # chromosome = svHelpers.clean_chromosomes([chromosome])[0]

        # pysam pileup does not allow pos1 = pos2. Thus we increase pos2
        # but will only consider the pileup at pos1.
        SAME_POS = False
        if position1 == position2:
            SAME_POS = True
            position2 += 1

        coverages = []
        for pileup in (bam_obj.pileup(chromosome,
                                      int(position1 - self.PYSAM_ZERO_BASED),
                                      int(position2 - self.PYSAM_ZERO_BASED))):
            # A: SAME_POS should now be handled by the break below
            #if ((pileup.pos + self.PYSAM_ZERO_BASED) >= position1) and ( ((pileup.pos + self.PYSAM_ZERO_BASED <= position2) and SAME_POS == False) or ((pileup.pos + self.PYSAM_ZERO_BASED < position2) and SAME_POS == True)):
            if ((pileup.pos + self.PYSAM_ZERO_BASED >= position1) and
                    (pileup.pos + self.PYSAM_ZERO_BASED <= position2)):

                coverage_counter = len(baseline_reads)
                if coverage_counter > 0:
                    for pileupread in pileup.pileups:
                        if pileupread.alignment.query_name in baseline_reads:
                            continue
                        else:
                            coverage_counter += 1
                else:
                    coverage_counter = (pileup.n)

                coverages.append(coverage_counter)

                # Only consider the pileup at pos1 if SAME_POS is True
                if (SAME_POS is True) and (pileup.pos + self.PYSAM_ZERO_BASED == position1):
                    coverages = [coverage_counter]
                    break

        cov_str = str(coverages)
        avg_cov = np.round(np.mean(coverages), 3)

        return (cov_str, avg_cov)
