import std.stdio;
import std.array: split;
import std.conv: to;
import std.algorithm: startsWith;
import std.getopt;

import dhtslib.vcf;
import dhtslib.bgzf;
import cadd;

enum Build{
    GRCh37,
    GRCh38
}
Build build = build.GRCh38;

void main(string[] args)
{
    auto getoptresult = getopt(args,
        "b|build", "Genome build used to generate CADD file (GRCh37 or GRCh38) default: GRCh38",&build,  
    );

    /// report options and exit
    if(getoptresult.helpWanted)
    {
        stderr.writefln("Usage: %s [options] <vcf(can be gzipped)> <cadd.tsv (can be gzippped)>", args[0]);
        defaultGetoptPrinter("cadd2vcf: Annotate vcf with CADD Report", getoptresult.options);
        return;
    }
    if(args.length < 3){
        stderr.writeln("Incorrect number of parameters");
        stderr.writefln("Usage: %s [options] <vcf(can be gzipped)> <cadd.tsv (can be gzippped)>", args[0]);
        defaultGetoptPrinter("cadd2vcf: Annotate vcf with CADD Report", getoptresult.options);
        return;
    }

    // open vcf readers and writer
    auto vcfreader = VCFReader(args[1]);
    auto vcfwriter = VCFWriter("-", vcfreader.vcfhdr);
    auto caddfile = BGZFile(args[2]);

    // add in new headers
    if(build == Build.GRCh37)
        addHeaderFields!true(vcfwriter);
    else
        addHeaderFields(vcfwriter);
    vcfwriter.writeHeader;

    // prime caddline
    auto caddEmpty = caddfile.empty;
    auto caddline = caddfile.front;

    // check for UCSC/chr prefix in vcf
    auto prefix = false;
    foreach (contig; vcfreader.vcfhdr.sequences)
    {
        if(contig.startsWith("chr")){
            prefix = true;
            break;
        }
    }

    // skip comment lines in caddfile
    while(caddline[0]=='#'){
        caddfile.popFront;
        caddEmpty = caddfile.empty;
        caddline = caddfile.front;
    }

    // loop
    foreach (VCFRecord rec; vcfreader)
    {
        rec.vcfheader = vcfwriter.getHeader;
        // if extra vcf records write and skip
        if(caddEmpty){
            vcfwriter.writeRecord(rec);
            continue;
        }


        auto caddfields = caddline.split("\t");
        // not extended record
        if(caddfields.length == 6){
            // process caddline and fix prefix in cadd record
            auto caddResult = processCADDLine!CADDAnno(caddfields);
            caddResult.pos--;
            if(prefix) caddResult.chr = "chr" ~ caddResult.chr; 
            if(
                    caddResult.chr != rec.chrom ||
                    caddResult.pos != rec.pos ||
                    caddResult.refAllele != rec.refAllele ||
                    caddResult.altAllele != rec.altAllelesAsArray[0]
                    ){
                vcfwriter.writeRecord(rec);
                continue;
            }
            // check that we are doing things correctly
            // TODO: Change to enforce for release checks?
            assert(caddResult.chr == rec.chrom);
            assert(caddResult.pos == rec.pos);
            assert(caddResult.refAllele == rec.refAllele);
            assert(caddResult.altAllele == rec.altAllelesAsArray[0]);

            // TODO: Change to enforce for release checks
            assert(rec.altAllelesAsArray.length == 1);

            addINFOFields(rec, caddResult);
            vcfwriter.writeRecord(rec);

        }else if(build == Build.GRCh37){
            // process caddline and fix prefix in cadd record
            auto caddResult = processCADDLine!CADDAnnoExtendedGRCH37(caddfields);
            caddResult.pos--;
            if(prefix) caddResult.chr = "chr" ~ caddResult.chr; 
            if(
                    caddResult.chr != rec.chrom ||
                    caddResult.pos != rec.pos ||
                    caddResult.refAllele != rec.refAllele ||
                    caddResult.altAllele != rec.altAllelesAsArray[0]
                    ){
                vcfwriter.writeRecord(rec);
                continue;
            }
            // check that we are doing things correctly
            // TODO: Change to enforce for release checks?
            assert(caddResult.chr == rec.chrom);
            assert(caddResult.pos == rec.pos);
            assert(caddResult.refAllele == rec.refAllele);
            assert(caddResult.altAllele == rec.altAllelesAsArray[0]);

            // TODO: Change to enforce for release checks
            assert(rec.altAllelesAsArray.length == 1);

            addINFOFields(rec, caddResult);
            vcfwriter.writeRecord(rec);
        }else{
            // process caddline and fix prefix in cadd record
            auto caddResult = processCADDLine!CADDAnnoExtended(caddfields);
            caddResult.pos--;
            if(prefix) caddResult.chr = "chr" ~ caddResult.chr; 
            if(
                    caddResult.chr != rec.chrom ||
                    caddResult.pos != rec.pos ||
                    caddResult.refAllele != rec.refAllele ||
                    caddResult.altAllele != rec.altAllelesAsArray[0]
                    ){
                vcfwriter.writeRecord(rec);
                continue;
            }
            // check that we are doing things correctly
            // TODO: Change to enforce for release checks?
            assert(caddResult.chr == rec.chrom);
            assert(caddResult.pos == rec.pos);
            assert(caddResult.refAllele == rec.refAllele);
            assert(caddResult.altAllele == rec.altAllelesAsArray[0]);

            // TODO: Change to enforce for release checks
            assert(rec.altAllelesAsArray.length == 1);

            addINFOFields(rec, caddResult);
            vcfwriter.writeRecord(rec);
        }

        // move caddfile range forward
        caddfile.popFront;
        caddEmpty = caddfile.empty;
        caddline = caddfile.front;
    }
}

