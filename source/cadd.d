module cadd;
import std.meta: AliasSeq, staticIndexOf;
import std.conv: to;
import dhtslib.vcf;


struct CADDAnno
{
	string chr;
	long pos;
	string refAllele;
	string altAllele;
	float rawScore;
	float phred;
}

struct CADDAnnoExtended
{
	string chr;
	long pos;
	string refAllele;
	string altAllele;
	string type;
	string length;
	string annoType;
	string Consequence;
    int consScore;
	string consDetail;
	float GC;
	float CpG;
	int motifECount;
	string motifEName;
	bool motifEHIPos;
	float motifEScoreChng;
	string oAA;
	string nAA;
	string geneID;
	string featureID;
	string geneName;
	string ccds;
	string intron;
	string exon;
	float cDNApos;
	float relcDNApos;
	float cDSpos;
	float relCDSpos;
	float protPos;
	float relProtPos;
	string domain;
	float dst2Splice;
	string dst2SplType;
	float minDistTSS;
	float minDistTSE;
	string SIFTcat;
	string SIFTval;
	float PolyPhenCat;
	float PolyPhenVal;
	float priPhCons;
	float mamPhCons;
	float verPhCons;
	float priPhyloP;
	float mamPhyloP;
	float verPhyloP;
	int bStatistic;
	int targetScan;
	float mirSVRScore;
	float mirSVRE;
	int mirSVRAln;
	float cHmmE1;
	float cHmmE2;
	float cHmmE3;
	float cHmmE4;
	float cHmmE5;
	float cHmmE6;
	float cHmmE7;
	float cHmmE8;
	float cHmmE9;
	float cHmmE10;
	float cHmmE11;
	float cHmmE12;
	float cHmmE13;
	float cHmmE14;
	float cHmmE15;
	float cHmmE16;
	float cHmmE17;
	float cHmmE18;
	float cHmmE19;
	float cHmmE20;
	float cHmmE21;
	float cHmmE22;
	float cHmmE23;
	float cHmmE24;
	float cHmmE25;
	float grepRS;
	float gerpRSpval;
	float grepN;
	float gerpS;
	float tOverlapMotifs;
	float motifDist;
	float EncodeH3K4me1Sum;
	float EncodeH3K4me1Max;
	float EncodeH3K4me2Sum;
	float EncodeH3K4me2Max;
	float EncodeH3K4me3Sum;
	float EncodeH3K4me3Max;
	float EncodeH3K9acMax;
	float EncodeH3K9acSum;
	float EncodeH3K9me3sum;
	float EncodeH3K9me3max;
	float EncodeH3K27acSum;
	float EncodeH3K27acMax;
    float EncodeH3K27me3Sum;
    float EncodeH3K27me3Max;
	float EncodeH3K36me3sum;
	float EncodeH3K36me3max;
	float EncodeH3K79me2sum;
	float EncodeH3K79me2max;
	float EncodeH4K20me1sum;
	float EncodeH4K20me1max;
	float EncodeH2AFZSum;
	float EncodeH2AFZMax;
	float EncodeDNaseSum;
	float EncodeDNaseMax;
	float EncodetotalRNASum;
	float EncodetotalRNAMax;
	float Grantham;
	float SpliceAIAccGain;
	float SpliceAIAccLoss;
	float SpliceAIDonGain;
	float SpliceAIDonLoss;
	float MMSpAcceptorIntron;
	float MMSpAcceptor;
	float MMSpExon;
	float MMSpDonor;
	float MMSpDonorIntron;
	float dist2Mutation;
	int freq100bp;
	int rare100bp;
	int sngl100bp;
	int freq1000bp;
	int rare1000bp;
	int sngl1000bp;
	int freq10000bp;
	int rare10000bp;
	int sngl10000bp;
	string EnsembleRegulatoryFeature;
	float dbscSNVDdaScore;
	float dbscSNVRfScore;
	int remapOverlapTF;
	int remapOverlapCL;
	float rawScore;
	float phred;
}

alias CADD_SIMPLE_FIELDS = 
AliasSeq!(
    "CADD_RawScore",
	"CADD_PHRED",
);

alias CADD_SIMPLE_DESCRIPTIONS = 
AliasSeq!(
    "Raw score from the model",
	"CADD PHRED Score",
);

alias CADD_EXTENDED_FIELDS = 
AliasSeq!(
	"CADD_Type",
	"CADD_Length",
	"CADD_AnnoType",
	"CADD_Consequence",
	"CADD_ConsScore ",
	"CADD_ConsDetail",
	"CADD_GC",
	"CADD_CpG",
	"CADD_motifECount ",
	"CADD_motifEName",
	"CADD_motifEHIPos",
	"CADD_motifEScoreChng",
	"CADD_oAA",
	"CADD_nAA",
	"CADD_GeneID",
	"CADD_FeatureID",
	"CADD_GeneName",
	"CADD_CCDS",
	"CADD_Intron",
	"CADD_Exon",
	"CADD_cDNApos",
	"CADD_relcDNApos",
	"CADD_CDSpos",
	"CADD_relCDSpos",
	"CADD_protPos",
	"CADD_relProtPos",
	"CADD_Domain",
	"CADD_Dst2Splice",
	"CADD_Dst2SplType",
	"CADD_minDistTSS",
	"CADD_minDistTSE",
	"CADD_SIFTcat",
	"CADD_SIFTval",
	"CADD_PolyPhenCat",
	"CADD_PolyPhenVal",
	"CADD_priPhCons",
	"CADD_mamPhCons",
	"CADD_verPhCons",
	"CADD_priPhyloP",
	"CADD_mamPhyloP",
	"CADD_verPhyloP",
	"CADD_bStatistic ",
	"CADD_targetScan ",
	"CADD_mirSVR-Score",
	"CADD_mirSVR-E",
	"CADD_mirSVR-Aln ",
	"CADD_cHmm_E1",
	"CADD_cHmm_E2",
	"CADD_cHmm_E3",
	"CADD_cHmm_E4",
	"CADD_cHmm_E5",
	"CADD_cHmm_E6",
	"CADD_cHmm_E7",
	"CADD_cHmm_E8",
	"CADD_cHmm_E9",
	"CADD_cHmm_E10",
	"CADD_cHmm_E11",
	"CADD_cHmm_E12",
	"CADD_cHmm_E13",
	"CADD_cHmm_E14",
	"CADD_cHmm_E15",
	"CADD_cHmm_E16",
	"CADD_cHmm_E17",
	"CADD_cHmm_E18",
	"CADD_cHmm_E19",
	"CADD_cHmm_E20",
	"CADD_cHmm_E21",
	"CADD_cHmm_E22",
	"CADD_cHmm_E23",
	"CADD_cHmm_E24",
	"CADD_cHmm_E25",
	"CADD_GerpRS",
	"CADD_GerpRSpval",
	"CADD_GerpN",
	"CADD_GerpS",
	"CADD_tOverlapMotifs",
	"CADD_motifDist",
	"CADD_EncodeH3K4me1-sum",
	"CADD_EncodeH3K4me1-max",
	"CADD_EncodeH3K4me2-sum",
	"CADD_EncodeH3K4me2-max",
	"CADD_EncodeH3K4me3-sum",
	"CADD_EncodeH3K4me3-max",
	"CADD_EncodeH3K9ac-sum",
	"CADD_EncodeH3K9ac-max",
	"CADD_EncodeH3K9me3-sum",
	"CADD_EncodeH3K9me3-max",
	"CADD_EncodeH3K27ac-sum",
	"CADD_EncodeH3K27ac-max",
	"CADD_EncodeH3K27me3-sum",
	"CADD_EncodeH3K27me3-max",
	"CADD_EncodeH3K36me3-sum",
	"CADD_EncodeH3K36me3-max",
	"CADD_EncodeH3K79me2-sum",
	"CADD_EncodeH3K79me2-max",
	"CADD_EncodeH4K20me1-sum",
	"CADD_EncodeH4K20me1-max",
	"CADD_EncodeH2AFZ-sum",
	"CADD_EncodeH2AFZ-max",
	"CADD_EncodeDNase-sum",
	"CADD_EncodeDNase-max",
	"CADD_EncodetotalRNA-sum",
	"CADD_EncodetotalRNA-max",
	"CADD_Grantham",
	"CADD_SpliceAI-acc-gain",
	"CADD_SpliceAI-acc-loss",
	"CADD_SpliceAI-don-gain",
	"CADD_SpliceAI-don-loss",
	"CADD_MMSp_acceptorIntron",
	"CADD_MMSp_acceptor",
	"CADD_MMSp_exon",
	"CADD_MMSp_donor",
	"CADD_MMSp_donorIntron",
	"CADD_Dist2Mutation",
	"CADD_Freq100bp",
	"CADD_Rare100bp",
	"CADD_Sngl100bp",
	"CADD_Freq1000bp",
	"CADD_Rare1000bp",
	"CADD_Sngl1000bp",
	"CADD_Freq10000bp",
	"CADD_Rare10000bp",
	"CADD_Sngl10000bp",
	"CADD_EnsembleRegulatoryFeature",
	"CADD_dbscSNV-ada_score",
	"CADD_dbscSNV-rf_score",
	"CADD_RemapOverlapTF",
	"CADD_RemapOverlapCL",
	"CADD_RawScore",
	"CADD_PHRED",
);

alias CADD_EXTENDED_DESCRIPTIONS =
AliasSeq!(
	"Event type (SNV, DEL, INS)",
	"Number of inserted/deleted bases",
	"CodingTranscript, Intergenic, MotifFeature, NonCodingTranscript, RegulatoryFeature, Transcript",
	"VEP consequence, priority selected by potential impact (default: UNKNOWN)",
	"Custom deleterious score assigned to Consequence",
	"Trimmed VEP consequence prior to simplification",
	"Percent GC in a window of +/- 75bp (default: 0.42)",
	"Percent CpG in a window of +/- 75bp (default: 0.02)",
	"Total number of overlapping motifs (default: 0)",
	"Name of sequence motif the position overlaps",
	"Is the position considered highly informative for an overlapping motif by VEP (default: 0)",
	"VEP score change for the overlapping motif site (default: 0)",
	"Reference amino acid (default: unknown)",
	"Amino acid of observed variant (default: unknown)",
	"ENSEMBL GeneID",
	"ENSEMBL feature ID (Transcript ID or regulatory feature ID)",
	"GeneName provided in ENSEMBL annotation",
	"Consensus Coding Sequence ID",
	"Intron number/Total number of exons",
	"Exon number/Total number of exons",
	"Base position from transcription start (default: 0*)",
	"Relative position in transcript (default: 0)",
	"Base position from coding start (default: 0*)",
	"Relative position in coding sequence (default: 0)",
	"Amino acid position from coding start (default: 0*)",
	"Relative position in protein codon (default: 0)",
	"Domain annotation inferred from VEP annotation (ncoils, sigp,lcompl, hmmpanther, ndomain = \"other named domain\") (default:UD)",
	"Distance to splice site in 20bp; positive: exonic, negative: intronic(default: 0)",
	"Closest splice site is ACCEPTOR or DONOR (default: unknown)",
	"Distance to closest Transcribed Sequence Start (TSS) (default:5.5)",
	"Distance to closest Transcribed Sequence End (TSE) (default: 5.5)",
	"SIFT category of change (default: UD)",
	"SIFT score (default: 0*)",
	"PolyPhen category of change (default: UD)",
	"PolyPhen score (default: 0*)",
	"Primate PhastCons conservation score (excl. human) (default: 0.0)",
	"Mammalian PhastCons conservation score (excl. human) (default: 0.0)",
	"Vertebrate PhastCons conservation score (excl. human) (default: 0.0)",
	"Primate PhyloP score (excl. human) (default: -0.029)",
	"Mammalian PhyloP score (excl. human) (default: -0.005)",
	"Vertebrate PhyloP score (excl. human) (default: 0.042)",
	"Background selection score (default: 800)",
	"targetscan (default: 0*)",
	"mirSVR-Score (default: 0*)",
	"mirSVR-E (default: 0)",
	"mirSVR-Aln (default: 0)",
	"Number of 48 cell types in chromHMM state E1_poised (default: 1.92*)",
	"Number of 48 cell types in chromHMM state E2_repressed(default: 1.92)",
	"Number of 48 cell types in chromHMM state E3_dead (default: 1.92)",
	"Number of 48 cell types in chromHMM state E4_dead (default: 1.92)",
	"Number of 48 cell types in chromHMM state E5_repressed(default: 1.92)",
	"Number of 48 cell types in chromHMM state E6_repressed(default: 1.92)",
	"Number of 48 cell types in chromHMM state E7_weak (default: 1.92)",
	"Number of 48 cell types in chromHMM state E8_gene (default: 1.92)",
	"Number of 48 cell types in chromHMM state E9_gene (default: 1.92)",
	"Number of 48 cell types in chromHMM state E10_gene (default: 1.92)",
	"Number of 48 cell types in chromHMM state E11_gene (default: 1.92)",
	"Number of 48 cell types in chromHMM state E12_distal (default: 1.92)",
	"Number of 48 cell types in chromHMM state E13_distal (default: 1.92)",
	"Number of 48 cell types in chromHMM state E14_distal (default: 1.92)",
	"Number of 48 cell types in chromHMM state E15_weak (default: 1.92)",
	"Number of 48 cell types in chromHMM state E16_tss (default: 1.92)",
	"Number of 48 cell types in chromHMM state E17_proximal(default: 1.92)",
	"Number of 48 cell types in chromHMM state E18_proximal(default: 1.92)",
	"Number of 48 cell types in chromHMM state E19_tss (default: 1.92)",
	"Number of 48 cell types in chromHMM state E20_poised (default: 1.92)",
	"Number of 48 cell types in chromHMM state E21_dead (default: 1.92)",
	"Number of 48 cell types in chromHMM state E22_repressed(default: 1.92)",
	"Number of 48 cell types in chromHMM state E23_weak (default: 1.92)",
	"Number of 48 cell types in chromHMM state E24_distal (default: 1.92)",
	"Number of 48 cell types in chromHMM state E25_distal (default: 1.92)",
	"Gerp element score (default: 0)",
	"Gerp element p-Value (default: 0)",
	"Neutral evolution score defined by GERP++ (default: 3.0)",
	"Rejected Substitution score defined by GERP++ (default: -0.2)",
	"Number of overlapping predicted TF motifs",
	"Reference minus alternate allele difference in nucleotide frequency within an predicted overlapping motif (default: 0)",
	"Sum of Encode H3K4me1 levels (from 13 cell lines) (default: 0.76)",
	"Maximum Encode H3K4me1 level (from 13 cell lines) (default: 0.37)",
	"Sum of Encode H3K4me2 levels (from 14 cell lines) (default: 0.73)",
	"Maximum Encode H3K4me2 level (from 14 cell lines) (default: 0.37)",
	"Sum of Encode H3K4me3 levels (from 14 cell lines) (default: 0.81)",
	"Maximum Encode H3K4me3 level (from 14 cell lines) (default: 0.38)",
	"Sum of Encode H3K9ac levels (from 13 cell lines) (default: 0.82)",
	"Maximum Encode H3K9ac level (from 13 cell lines) (default: 0.41)",
	"Sum of Encode H3K9me3 levels (from 14 cell lines) (default: 0.81)",
	"Maximum Encode H3K9me3 level (from 14 cell lines) (default: 0.38)",
	"Sum of Encode H3K27ac levels (from 14 cell lines) (default: 0.74)",
	"Maximum Encode H3K27ac level (from 14 cell lines) (default: 0.36)",
	"Sum of Encode H3K27me3 levels (from 14 cell lines) (default: 0.93)",
	"Maximum Encode H3K27me3 level (from 14 cell lines) (default: 0.47)",
	"Sum of Encode H3K36me3 levels (from 10 cell lines) (default: 0.71)",
	"Maximum Encode H3K36me3 level (from 10 cell lines) (default: 0.39)",
	"Sum of Encode H3K79me2 levels (from 13 cell lines) (default: 0.64)",
	"Maximum Encode H3K79me2 level (from 13 cell lines) (default: 0.34)",
	"Sum of Encode H4K20me1 levels (from 11 cell lines) (default: 0.88)",
	"Maximum Encode H4K20me1 level (from 11 cell lines) (default: 0.47)",
	"Sum of Encode H2AFZ levels (from 13 cell lines) (default: 0.9)",
	"Maximum Encode H2AFZ level (from 13 cell lines) (default: 0.42)",
	"Sum of Encode DNase-seq levels (from 12 cell lines) (default: 0.0)",
	"Maximum Encode DNase-seq level (from 12 cell lines) (default: 0.0)",
	"Sum of Encode totalRNA-seq levels (from 10 cell lines always minus and plus strand) (default: 0.0)",
	"Maximum Encode totalRNA-seq level (from 10 cell lines, minus and plus strand separately) (default: 0.0)",
	"Grantham score: oAA,nAA (default: 0*)",
	"Masked SpliceAI acceptor gain score (default: 0*)",
	"Masked SpliceAI acceptor loss score (default: 0)",
	"Masked SpliceAI donor gain score (default: 0)",
	"Masked SpliceAI donor loss score (default: 0)",
	"MMSplice acceptor intron (intron 3’) score (default: 0)",
	"MMSplice acceptor score (default: 0)",
	"MMSplice exon score (default: 0)",
	"MMSplice donor score (default: 0)",
	"MMSplice donor intron (intron 5’) )score (default: 0)",
	"Distance between the closest BRAVO SNV up and downstream (position itself excluded) (default: 0*)",
	"Number of frequent (MAF > 0.05) BRAVO SNV in 100 bp window nearby (default: 0)",
	"Number of rare (MAF < 0.05) BRAVO SNV in 100 bp window nearby (default: 0)",
	"Number of single occurrence BRAVO SNV in 100 bp window nearby (default: 0)",
	"Number of frequent (MAF > 0.05) BRAVO SNV in 1000 bp window nearby (default: 0)",
	"Number of rare (MAF < 0.05) BRAVO SNV in 1000 bp window nearby (default: 0)",
	"Number of single occurrence BRAVO SNV in 1000 bp window nearby (default: 0)",
	"Number of frequent (MAF > 0.05) BRAVO SNV in 10000 bp window nearby (default: 0)",
	"Number of rare (MAF < 0.05) BRAVO SNV in 10000 bp window nearby (default: 0)",
	"Number of single occurrence BRAVO SNV in 10000 bp window nearby (default: 0)",
	"Matches in the Ensemble Regulatory Built (similar to annotype) (default: NA)",
	"Adaboost classifier score from dbscSNV (default: 0*)",
	"Random forest classifier score from dbscSNV (default: 0*)",
	"Remap number of different transcription factors binding (default: -0.5)",
	"Remap number of different transcription factor - cell line combinations binding (default: -0.5)",
	"Raw score from the model",
	"CADD PHRED Score",
);

alias VcfTypes = AliasSeq!(float, int, bool, string);
enum VcfTypeIndex(T) = staticIndexOf!(T, VcfTypes);
/// See https://samtools.github.io/hts-specs/SAMv1.pdf sec 1.5
alias TypeStrings = AliasSeq!("Float","Integer","Flag","String");

auto processCADDLine(T)(string[] fields){
    alias memberTypes = AliasSeq!(typeof(T.tupleof));
    alias memberNames = AliasSeq!(T.tupleof);
	assert(fields.length == memberTypes.length);
    assert(fields.length == memberNames.length);
    T ret;
    static foreach(i,name;memberNames){
        mixin("ret."~name.stringof~"=fields[" ~ i.to!string ~ "].to!"~memberTypes[i].stringof~";");
    }
	return ret;
}

void addHeaderFields(ref VCFWriter writer){
    alias memberTypes = AliasSeq!(typeof(CADDAnnoExtended.tupleof));
    static foreach (i,m; memberTypes[4..$])
    {
        mixin("writer.addTag!\"INFO\"(" ~ CADD_EXTENDED_FIELDS[i].stringof ~ ",\"A\"," ~ TypeStrings[VcfTypeIndex!m].stringof ~ "," ~ CADD_EXTENDED_DESCRIPTIONS[i].stringof ~");");    
    }
}

void addINFOFields(T)(ref VCFRecord rec, ref T anno){
    alias memberNames = AliasSeq!(T.tupleof);
    static foreach (i,name; memberNames[4..$])
    {
        static if(is(T == CADDAnno)){
            mixin("rec.addInfo(" ~ CADD_SIMPLE_FIELDS[i].stringof ~ ",anno." ~ name.stringof ~ ");");
        }else{
            mixin("rec.addInfo(" ~ CADD_EXTENDED_FIELDS[i].stringof ~ ",anno." ~ name.stringof ~ ");");
        }
    }
}

// string processCADDLineMixin(T,types,string var)(){
//     string ret = "return " ~ T.stringof ~"(";
//     static foreach (i,name; names)
//     {
//         ret ~= var~
//     }
// }