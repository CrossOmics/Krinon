#include "Parameters.h"

const int SJDB_PADDING_LENGTH = 20;

/**
 * TODO: A bunch of values are hard-coded (I cannot comment if they have been inherited from STAR,
 *       or tuned by some other logic), and they make reasoning about things very difficult.
 *       The names also overlap (there are multiple `kMerSize` values for example).
 *       I propose we bring them all together into a config file and set them as default values
 *       here so we can get a look at them from just one place.
 *       My current file is by no means complete, so we'll have to sweep through the code at some
 *       point before it gets out of hand.
 */

namespace rna {
    Parameters::Parameters() : 
        program("rnaAligner")
    {
        // First, discover where we are even running
        char result[PATH_MAX];
        if (readlink("/proc/self/exe", result, PATH_MAX) > 0) {
            m_binaryDir = std::filesystem::path(result).parent_path();
            PLOG_INFO << "Binary running on path: " << m_binaryDir;
        }
        else {
            PLOG_FATAL << "Couldn't discover running binary location";
            exit(-1);
        }

        // Load the config file
        auto configPath = m_binaryDir / std::filesystem::path("config.yaml");
        if (std::filesystem::exists(configPath)) {
            PLOG_INFO << "Found configuration file at: " << configPath;
            m_globalConfig = YAML::LoadFile(configPath.c_str());
        }
        else {
            PLOG_FATAL << "Could not locate configuration file!";
            exit(-1);
        }
        auto GENOME_INDEX_CONFIG = m_globalConfig["genome"]["GenomeIndexConfig"];
        auto GENOME_PREFIX_CONFIG = m_globalConfig["genome"]["GenomeIndexPrefixConfig"];
        auto READ_ALIGNER_CONFIG = m_globalConfig["readAlign"]["ReadAligner"];
        auto SEED_MAPPING_CONFIG = m_globalConfig["readAlign"]["SeedMappingConfig"];
        auto STITCHING_CONFIG = m_globalConfig["readAlign"]["StitchingConfig"];

        program.add_argument("--mode")
            .help("the mode to run: GenomeGenerate, ReadAlign, GenomeGenerateAndReadAlign, test")
            .default_value(std::string(""));
        program.add_argument("--genomeBinSize")
            .help("log2 of the size of each bin in the genome index, default 18")
            .default_value(18)
            .scan<'i',int>();
        program.add_argument("--genomeFile")
            .help("the reference genome file in fasta format")
            .default_value("");
        program.add_argument("--readFile")
            .help("the read file in fastq format")
            .default_value("");
        program.add_argument("--readFile2")
            .help("the second read file for paired-end reads, optional")
            .default_value("");
        program.add_argument("--readType")
            .help("the type of reads: single or paired, default single")
            .default_value(std::string("single"));
        program.add_argument("--gtfFile")
            .help("the gene annotation file in gtf format, optional")
            .default_value("");
        program.add_argument("--threads")
            .help("number of threads to use , default 1")
            .default_value(1)
            .scan<'i',int>();
        program.add_argument("--outPutDir")
            .help("the directory to save the output files, default current directory")
            .default_value(std::string("./"));
        program.add_argument("--genomeGenerateFileStoreDir")
            .help("the directory to save the genome index files, default current directory")
            .default_value(std::string("./"));
        
        auto &alignerGroup = program.add_group("Aligner Parameters");

        alignerGroup.add_argument("--inputBufferSize")
            .help("Number of reads to load from input file as a single chunk")
            .default_value(READ_ALIGNER_CONFIG["inputBufferSize"].as<int>())
            .scan<'i', int>();
        alignerGroup.add_argument("--outputBufferSize")
            .help("Number of characters per each output chunk for writing")
            .default_value(READ_ALIGNER_CONFIG["outputBufferSize"].as<int>())
            .scan<'i', int>();
        
        auto &sjdbGroup = program.add_group("Splice Junction DB Parameters");

        sjdbGroup.add_argument("--kMerSize")
        // TODO: This value is later changed to 10, as noted in the config file. Which on is it?
            .help("k-mer size used in the genome index, default 14")
            .default_value(14)
            .scan<'i',int>();
        sjdbGroup.add_argument("--sjdbOverhang")
            .help("the length of the donor/acceptor sequence around each splice junction")
            .default_value(GENOME_INDEX_CONFIG["sjdbOverhang"].as<int>())
            .scan<'i',int>();
        sjdbGroup.add_argument("--limitSjdbInsertN")
            .help("the maximum number of splice junctions to insert into the genome index")
            .default_value(GENOME_INDEX_CONFIG["limitSjdbInsertN"].as<int>())
            .scan<'i',int>();
        sjdbGroup.add_argument("--sjdbGTFfeatureExon")
            .help("the feature type in the GTF file to use as exons, default exon")
            .default_value(std::string("exon"));
        sjdbGroup.add_argument("--sjdbGTFTagExonParentTranscriptId")
            .help("the attribute type in the GTF file to use as transcript ID, default transcript_id")
            .default_value(std::string("transcript_id"));
        sjdbGroup.add_argument("--sjdbGTFChrPrefix")
            .help("the prefix of chromosome names in the GTF file, default empty")
            .default_value(std::string(""));
        sjdbGroup.add_argument("--sjdbGTFTagExonParentGene")
            .help("the attribute type in the GTF file to use as parent gene ID, default gene_id")
            .default_value(std::string("gene_id"));
        sjdbGroup.add_argument("--sjdbGTFTagExonParentGeneName")
            .nargs(argparse::nargs_pattern::any)
            .help("the attribute type in the GTF file to use as parent gene name, default gene_name")
            .default_value(std::vector<std::string>{std::string("gene_name")});
        sjdbGroup.add_argument("--sjdbGTFTagExonParentGeneType")
            .nargs(argparse::nargs_pattern::any)
            .help("the attribute type in the GTF file to use as parent gene type, default gene_type")
            .default_value(std::vector<std::string>{std::string("gene_type"),std::string ("gene_biotype")});
        
        auto &seedMappingGroup = program.add_group("Seed Mapping Parameters");

        seedMappingGroup.add_argument("--minSplitLength")
            .help("the minimum length of each segment when splitting the read for alignment")
            .default_value(SEED_MAPPING_CONFIG["minSplitLength"].as<int>())
            .scan<'i',int>();
        seedMappingGroup.add_argument("--maxSeedPerRead")
            .help("the maximum number of seeds for each read")
            .default_value(SEED_MAPPING_CONFIG["maxSeedPerRead"].as<int>())
            .scan<'i',int>();
        
        auto &stitchingGroup = program.add_group("Stitching Parameters");

        /**
         * TODO: Some parameter names are not he same in the config file (and thus the code),
         *       e.g. `maxWindowNum` becomes just `windowNum`.
         */

        stitchingGroup.add_argument("--maxAnchorRep")
            .help("the maximum number of times an anchor can appear in the genome")
            .default_value(STITCHING_CONFIG["maxAnchorRep"].as<int>())
            .scan<'i',int>();
        stitchingGroup.add_argument("--winBinSize")
            .help("log2 of the size of each window for stitching")
            .default_value(STITCHING_CONFIG["winBinSizeLog"].as<int>())
            .scan<'i',int>();
        stitchingGroup.add_argument("--winAnchorDistBins")
            .help("the number of bins for the distance between anchors in adjacent windows")
            .default_value(STITCHING_CONFIG["winAnchorDistBins"].as<int>())
            .scan<'i',int>();
        stitchingGroup.add_argument("--winFlankSize")
            .help("the size of the flanking region for each window")
            .default_value(STITCHING_CONFIG["flankSize"].as<int>())
            .scan<'i',int>();
        stitchingGroup.add_argument("--maxWindowNum")
            .help("the maximum number of windows to use for stitching")
            .default_value(STITCHING_CONFIG["maxWindows"].as<int>())
            .scan<'i',int>();
        stitchingGroup.add_argument("--maxSeedPerWindow")
            .help("the maximum number of seeds in each window")
            .default_value(STITCHING_CONFIG["maxSeedPerWindows"].as<int>())
            .scan<'i',int>();
        stitchingGroup.add_argument("--maxAlignmentRep")
            .help("the maximum number of times an alignment can appear in the genome")
            .default_value(STITCHING_CONFIG["maxRep"].as<int>())
            .scan<'i',int>();
        stitchingGroup.add_argument("--maxExonNum")
            .help("the maximum number of exons in a transcript")
            .default_value(STITCHING_CONFIG["maxExons"].as<int>())
            .scan<'i',int>();
        stitchingGroup.add_argument("--maxTranscriptStoredNum")
            .help("the maximum number of transcripts to store during stitching")
            .default_value(STITCHING_CONFIG["transcriptStoredMax"].as<int>())
            .scan<'i',int>();
        stitchingGroup.add_argument("--outFilterMultimapMax")
            .help("the maximum number of multiple alignments allowed for each read")
            .default_value(STITCHING_CONFIG["outFilterMultimapMax"].as<int>())
            .scan<'i',int>();
        stitchingGroup.add_argument("--maxMismatch")
            .help("the maximum number of mismatches allowed in each read")
            .default_value(STITCHING_CONFIG["maxMismatch"].as<int>())
            .scan<'i',int>();
        /**
         * TODO: This scoring should be a `double` value, right?
         */
        stitchingGroup.add_argument("--multimapScoreRange")
            .help("the score range for multiple alignments")
            .default_value(STITCHING_CONFIG["multimapScoreRange"].as<int>())
            .scan<'i',int>();
        stitchingGroup.add_argument("--outFilterScoreMinOverLread")
            .help("the minimum ratio of alignment score to read length")
            .default_value(STITCHING_CONFIG["outFilterScoreMinOverLRead"].as<double>())
            .scan<'g',double>();
        stitchingGroup.add_argument("--outFilterMatchMinOverLread")
            .help("the minimum ratio of matched bases to read length")
            .default_value(STITCHING_CONFIG["outFilterMatchMinOverLRead"].as<double>())
            .scan<'g',double>();
        
        /**
         * TODO: I am extremely confused by the values here, not what they are/mean, but where
         * they are actually finally settled on.
         * Here they get a new default value, later in `StitchingManagement.h` they are given
         * default values multiple times, and then forced to be constant (I don't if even setting
         * them here actually changes something).
         */

        auto &scoringGroup = program.add_group("Scoring Parameters");

        //scoring
        scoringGroup.add_argument("--MATCH_SCORE")
            .help("the score for a match, default 1")
            .default_value(1)
            .scan<'i',int>();
        scoringGroup.add_argument("--MISMATCH_PENALTY")
            .help("the penalty for a mismatch, default -1")
            .default_value(-1)
            .scan<'i',int>();
        scoringGroup.add_argument("--GAP_OPEN_PENALTY")
            .help("the penalty for opening a gap, default 0")
            .default_value(0)
            .scan<'i',int>();
        scoringGroup.add_argument("--DEL_OPEN_PENALTY")
            .help("the penalty for opening a deletion, default -2")
            .default_value(-2)
            .scan<'i',int>();
        scoringGroup.add_argument("--DEL_EXTEND_PENALTY")
            .help("the penalty for extending a deletion, default -2")
            .default_value(-2)
            .scan<'i',int>();
        scoringGroup.add_argument("--INS_OPEN_PENALTY")
            .help("the penalty for opening an insertion, default -2")
            .default_value(-2)
            .scan<'i',int>();
        scoringGroup.add_argument("--INS_EXTEND_PENALTY")
            .help("the penalty for extending an insertion, default -2")
            .default_value(-2)
            .scan<'i',int>();
        scoringGroup.add_argument("--SCORE_STITCH_SJ_SHIFT")
            .help("the maximum score reduction while searching for splice junction boundaries, default 1")
            .default_value(1)
            .scan<'i',int>();
        scoringGroup.add_argument("--SCORE_GAP_GCAG")
            .help("the penalty for a GC-AG or CT-GC splice junction, default -4")
            .default_value(-4)
            .scan<'i',int>();
        scoringGroup.add_argument("--SCORE_GAP_ATAC")
            .help("the penalty for an AT-AC or GT-AT splice junction, default -8")
            .default_value(-8)
            .scan<'i',int>();
        scoringGroup.add_argument("--SCORE_GAP_NON_CANONICAL")
            .help("the penalty for a non-canonical splice junction, default -8")
            .default_value(-8)
            .scan<'i',int>();
        scoringGroup.add_argument("--SCORE_ANNOTATED_SJ")
            .help("the bonus score for an annotated splice junction, default 2")
            .default_value(2)
            .scan<'i',int>();
        scoringGroup.add_argument("--MAX_SJ_REPEAT_SEARCH")
            .help("the maximum number of bases to search for repeats around a splice junction, default 255")
            .default_value(255)
            .scan<'i',int>();
        scoringGroup.add_argument("--MIN_INTRON_LENGTH")
            .help("the minimum length of an intron, default 21")
            .default_value(21)
            .scan<'i',int>();
        scoringGroup.add_argument("--MAX_INTRON_LENGTH")
            .help("the maximum length of an intron, default 2147483647 (no limits)")
            .default_value(2147483647)
            .scan<'i',int>();
        scoringGroup.add_argument("--MAX_MISMATCH_FOR_SJ")
            .nargs(4)
            .help("the maximum number of mismatches allowed for different splice junction types: non-canonical, GT-AG, GC-AG, AT-AC. -1 means no limit. Default 0 -1 0 0")
            .default_value(std::vector<int>{0, -1, 0, 0})
            .scan<'i',int>();
    }

    int Parameters::process(int argc, char* argv[]){
        try {
            program.parse_args(argc, argv);
        }
        catch (const std::exception& err) {
            std::cerr << err.what() << std::endl;
            std::cerr << program;
            return 1;
        };
        mode = program.get<std::string>("--mode");
        genomeBinSize = program.get<int>("--genomeBinSize");
        genomeFile = program.get<std::string>("--genomeFile");
        readFile = program.get<std::string>("--readFile");
        readFile2 = program.get<std::string>("--readFile2");
        isPaired = program.get<std::string>("--readType") == "paired";
        gtfFile = program.get<std::string>("--gtfFile");
        threads = program.get<int>("--threads");
        outPutDir = program.get<std::string>("--outPutDir");
        genomeGenerateFileStoreDir = program.get<std::string>("--genomeGenerateFileStoreDir");
        outLogFile = std::ofstream(outPutDir + "/log.txt");
        int sjdbOverhang = program.get<int>("--sjdbOverhang");
        int sjdbLength = sjdbOverhang * 2 + SJDB_PADDING_LENGTH;
        genomeIndexConfig = GenomeIndexConfig{
            .kMerSize = program.get<int>("--kMerSize"),
            .insertSJ = !gtfFile.empty(),
            .sjdbOverhang = sjdbOverhang,
            .sjdbLength = sjdbLength,
            .limitSjdbInsertN =program.get<int>("--limitSjdbInsertN"),
        };
        gtfConfig = GTFConfig{
            .sjdbGTFfeatureExon = program.get<std::string>("--sjdbGTFfeatureExon"),
            .sjdbGTFTagExonParentTranscriptId = program.get<std::string>("--sjdbGTFTagExonParentTranscriptId"),
            .sjdbGTFChrPrefix = program.get<std::string>("--sjdbGTFChrPrefix"),
            .sjdbGTFTagExonParentGene = program.get<std::string>("--sjdbGTFTagExonParentGene"),
            .sjdbGTFTagExonParentGeneName = program.get<std::vector<std::string>>("--sjdbGTFTagExonParentGeneName"),
            .sjdbGTFTagExonParentGeneType = program.get<std::vector<std::string>>("--sjdbGTFTagExonParentGeneType"),
            .sjdbOverhang = sjdbOverhang,
            .sjdbLength = sjdbLength,
            .limitSjdbInsertN =program.get<int>("--limitSjdbInsertN"),
        };
        seedMappingConfig = SeedMappingConfig{
            .minSplitLength = program.get<int>("--minSplitLength"),
            .maxSeedPerRead = program.get<int>("--maxSeedPerRead"),
        };
        stitchConfig = StitchingConfig{
            .maxAnchorRep = program.get<int>("--maxAnchorRep"),
            .winBinSizeLog = program.get<int>("--winBinSize"),
            .winAnchorDistBins = program.get<int>("--winAnchorDistBins"),
            .flankSize = program.get<int>("--winFlankSize"),
            .maxWindows = program.get<int>("--maxWindowNum"),
            .maxSeedPerWindows = program.get<int>("--maxSeedPerWindow"),
            .maxRep = program.get<int>("--maxAlignmentRep"),
            .maxExons = program.get<int>("--maxExonNum"),
            .transcriptStoredMax = program.get<int>("--maxTranscriptStoredNum"),
            .outFilterMultimapMax = program.get<int>("--outFilterMultimapMax"),
            .maxMismatch = program.get<int>("--maxMismatch"),
            .multimapScoreRange = program.get<int>("--multimapScoreRange"),
            .outFilterScoreMinOverLRead = program.get<double>("--outFilterScoreMinOverLread"),
            .outFilterMatchMinOverLRead = program.get<double>("--outFilterMatchMinOverLread"),
        };
        auto maxMismatchForSJVec = program.get<std::vector<int>>("--MAX_MISMATCH_FOR_SJ");
        stitchingScoreConfig = StitchingScoreConfig{
            .MATCH_SCORE = program.get<int>("--MATCH_SCORE"),
            .MISMATCH_PENALTY = program.get<int>("--MISMATCH_PENALTY"),
            .GAP_OPEN_PENALTY = program.get<int>("--GAP_OPEN_PENALTY"),
            .DEL_OPEN_PENALTY = program.get<int>("--DEL_OPEN_PENALTY"),
            .DEL_EXTEND_PENALTY = program.get<int>("--DEL_EXTEND_PENALTY"),
            .INS_OPEN_PENALTY = program.get<int>("--INS_OPEN_PENALTY"),
            .INS_EXTEND_PENALTY = program.get<int>("--INS_EXTEND_PENALTY"),
            .SCORE_STITCH_SJ_SHIFT = program.get<int>("--SCORE_STITCH_SJ_SHIFT"),
            .SCORE_GAP_GCAG = program.get<int>("--SCORE_GAP_GCAG"),
            .SCORE_GAP_ATAC = program.get<int>("--SCORE_GAP_ATAC"),
            .SCORE_GAP_NON_CANONICAL = program.get<int>("--SCORE_GAP_NON_CANONICAL"),
            .SCORE_ANNOTATED_SJ = program.get<int>("--SCORE_ANNOTATED_SJ"),
            .MAX_SJ_REPEAT_SEARCH = program.get<int>("--MAX_SJ_REPEAT_SEARCH"),
            .MIN_INTRON_LENGTH = program.get<int>("--MIN_INTRON_LENGTH"),
            .MAX_INTRON_LENGTH = program.get<int>("--MAX_INTRON_LENGTH"),
            .MAX_MISMATCH_FOR_SJ = {maxMismatchForSJVec[0], maxMismatchForSJVec[1], maxMismatchForSJVec[2], maxMismatchForSJVec[3]},
        };

        return 0;
    }
}