// ==========================================================================
//                           Mapping SMRT reads 
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: cxpan <chenxu.pan@fu-berlin.de>
// ==========================================================================

#include <csignal>
#include "mapper.h"
#include "dp_verifier.h"

using namespace seqan;


    seqan::ArgumentParser::ParseResult
    parseCommandLine(Options & options, int argc, char const ** argv)
    {
        // Setup ArgumentParser.
        seqan::ArgumentParser parser("pacMapper");
        // Set short description, version, and date.
        setShortDescription(parser, "Alignment of SMRT sequencing read");
        setVersion(parser, "1.0");
        setDate(parser, "May 2017");

        // Define usage line and long description.
        addUsageLine(parser,
                     "[\\fIOPTIONS\\fP] \"\\fIread.fa\\fP\" \"\\fIgnome.fa\\fP\"");
        addDescription(parser,
                       "Program for mapping raw SMRT sequencing reads to reference genome.");

        // Argument.
        addArgument(parser, seqan::ArgParseArgument(
            seqan::ArgParseArgument::INPUT_FILE, "read"));
        setHelpText(parser, 0, "Reads file .fa, .fasta");

        addArgument(parser, seqan::ArgParseArgument(
            seqan::ArgParseArgument::INPUT_FILE, "genome", true));
        setHelpText(parser, 1, "Reference file .fa, .fasta");

        addSection(parser, "Mapping Options");
        addOption(parser, seqan::ArgParseOption(
            "o", "output", "choose output file.",
             seqan::ArgParseArgument::STRING, "STR"));
        addOption(parser, seqan::ArgParseOption(
            "s", "sensitive", "Sensitive mode. Default closed",
             seqan::ArgParseArgument::STRING, "STR"));

        // Add Examples Section.
        addTextSection(parser, "Examples");
        addListItem(parser,
                    "\\fBpacMapper\\fP \\fB-U\\fP \\fIchr1.fa reads.fa\\fP",
                    "Print version of \"rd\"");
        addListItem(parser,
                    "\\fBpacMapper\\fP \\fB-L\\fP \\fB-i\\fP \\fI3\\fP "
                    "\\fIchr1.fa reads.fa\\fP",
                    "Print \"\" with every third character "
                    "converted to upper case.");

        // Parse command line.
        seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

        if (res != seqan::ArgumentParser::PARSE_OK)
            return res;

        getOptionValue(options.oPath, parser, "output");

        seqan::getArgumentValue(options.rPath, parser, 0);
        seqan::getArgumentValue(options.gPath, parser, 1);


        return seqan::ArgumentParser::PARSE_OK;

    }

    
int main(int argc, char const ** argv)
{
    std::cerr << "Encapsulated version: Mapping reads efficiently" << std::endl;
    (void)argc;
    // Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    double t=sysTime();
    Mapper<> mapper(options);

    double time = sysTime();
    mapper.createIndex();
    std::cerr << "done1\n";
    resize(mapper.hits(), length(mapper.reads()));
    resize(mapper.cords(), length(mapper.reads()));
    std::cerr << "done2\n";

    rawMapAll<typename MapperBase<>::DefaultAlphabet, typename MapperBase<>::DefaultShape>(mapper.index(), mapper.reads(), mapper.genomes(),
                         _DefaultMapParm, mapper.hits(), mapper.cords());
    auto read_map_coordinates = mapper.coordsAsVector();

    int grid_size = 1000;

    Score<int, Simple> scoreScheme(4, -5, -1, -11);

    for (int read = 0; read < length(read_map_coordinates); ++read)
    {
        std::cout << "\n\n#### Read " << read << "\n";
        for (int mapping = 0; mapping < length(read_map_coordinates[read]); ++mapping)
        {
            std::cout << "# Mapping " << mapping << "\n";
            PacBioBandGenerator bandGenerator{};
            computeGrid(bandGenerator, read_map_coordinates[read][mapping], grid_size);

            // We got the blocks here.
            // Now submit an alignment with the block generation as function?
            DPSettings<decltype(scoreScheme), DPVerifierTraits> dpSettings{scoreScheme};
            dpSettings.bandScheme = bandGenerator;
        }
    }


    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;

    return 0;
}
