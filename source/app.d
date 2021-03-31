import std.stdio;
import std.array: split;
import std.conv: to;
import std.algorithm: startsWith;

import dhtslib.vcf;
import dhtslib.bgzf;
import cadd;

void main(string[] args)
{
	// open vcf readers and writer
	auto vcfreader = VCFReader(args[1]);
	auto vcfwriter = VCFWriter("-", vcfreader.vcfhdr);
	auto caddfile = BGZFile(args[2]);

	// add in new headers
	addHeaderFields(vcfwriter);
	vcfwriter.writeHeader;

	// prime caddline
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
		caddline = caddfile.front;
	}

	// loop
	foreach (VCFRecord rec; vcfreader)
	{
		rec.vcfheader = vcfwriter.getHeader;
		// if extra vcf records write and skip
		if(caddfile.empty){
			vcfwriter.writeRecord(rec);
			continue;
		}


		auto caddfields = caddline.split("\t");
		// not extended record
		if(caddfields.length == 6){
			// process caddline and fix prefix in cadd record
			auto caddResult = processCADDLine!CADDAnno(caddfields);
			if(prefix) caddResult.chr = "chr" ~ caddResult.chr; 

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

		}else if(caddfields.length == 134){
			// process caddline and fix prefix in cadd record
			auto caddResult = processCADDLine!CADDAnnoExtended(caddfields);
			if(prefix) caddResult.chr = "chr" ~ caddResult.chr; 

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
			throw new Exception("malformed cadd line: " ~ caddline);
		}

		// move caddfile range forward
		caddfile.popFront;
		caddline = caddfile.front;
	}
}

