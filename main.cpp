#include "main.hpp"
#include "kmerUtilities.hpp"
#include "smithWaterman.hpp"
#include "hash.hpp"
#include "serialiseKmersMap.hpp"
#include <map>
#include <cstring>
#include <filesystem>
#include <iostream>

using std::cout;
using std::endl;
using std::cerr;
std::vector<std::string> sampleSequences;
std::vector<std::tuple<int, int, int, int>> results;  // Vector to store the matching results
std::mutex resultsMutex; // Mutex to protect concurrent writes to the results vector
std::vector<std::string> fileNames;
std::vector<std::tuple<uint32_t, uint32_t, uint32_t, double, std::string>> resultsSW; // match results after Smith Waterman alignment

// some global default values (overide with cmd line arguments)
unsigned int KMERSIZE; // seed kmer size
unsigned int MIN_MATCHES;   // minimum number of seed matches required
unsigned int SEED_SKIP; // seed skipping distancee
unsigned int THREADS; // number of trheads (default is HW concurrency)
int CUTOFF; // Smith Waterman cutoff (minimum number of aligned bases)
std::string REFERENCE; // file name for reference file to match against
std::string READS; // file name for reads to be matched by KISS
std::string OUTPUT_DIR; // directory name for results
uint64_t kmer_encoding_table[256]; // for packing a kmer into uint64
bool MATCH_ALL; // find all mathces flag
bool GET_SW; // compute Smith Waterman flag


/*
  GET BASE FILE NAME
  Extract the base filename without the path and extension
 */
std::string getBaseFilename(const std::string& filePath) {
    std::size_t lastSlash = filePath.find_last_of("/\\");
    std::size_t lastDot = filePath.find_last_of(".");
    if (lastDot != std::string::npos && (lastSlash == std::string::npos || lastDot > lastSlash)) {
        return filePath.substr(lastSlash + 1, lastDot - lastSlash - 1);
    }
    return filePath.substr(lastSlash + 1);
}



/*
  WRITE MATCHES TO FILE
 This will have to be replaced with a writer thread if we split the input file (being filetered) into two output files, one with the matching reads, and one with the clear reads.  So every read is filtered and immediately written out as required by the program parameters.
 */
void writeMatchesToFile(const std::string& filename,
                const std::vector<std::tuple<uint32_t, uint32_t, uint32_t, double, std::string>>& resultsSW,
                QueryTable ref_ID_Table) {
   
    // Extract the base filename
    std::string baseFilename = getBaseFilename(filename);

    // Create a new filename with the ".out" extension
    std::string newFilename = OUTPUT_DIR + "/" + baseFilename + ".out";

    // Open the new file for writing
    std::ofstream outFile(newFilename);

    // Check if the file was successfully opened
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << newFilename << std::endl;
        return;
    }
     //    ref_ID_Table.printKeys();
        
    // Print the resultsSW vector
    int n = 0;
    int unique=0;
    int prev=-1;
    uint32_t k;
    for (const auto& result : resultsSW) {
        uint32_t sampleIndex, refBegin, refEnd;
        double swScore;
        std::string cigar;
        std::tie(sampleIndex, refBegin, refEnd, swScore, cigar) = result;
        if (sampleIndex!=prev) {
            unique++;
            prev = sampleIndex;
        }
       
        // key is reference# in fasta file
        // note: the queryTable struct is define in serialiseKmersMap (for some unkown reason... :-) )
        uint32_t key = ref_ID_Table.query(refBegin);
        if (key == -1) {
            std::cout << "refrenceIndex = " << refBegin << " Ref_ID not found." << std::endl;
        }

        k = ref_ID_Table.getKeyFromIndex(key);
        std::ostringstream oss;
        oss << std::setw(4) << std::right << ++n << " "                  // match number
            << std::setw(9) << std::right << (sampleIndex+1) << " "      // read number, 1-based
            << std::setw(5) << std::right << key+1 << " "                  // ref sequence# 1-based
            << std::setw(9) << std::right << refBegin << " "   // match start position in ref sequence
            << std::setw(9) << std::right << refEnd << " "   //   match end position in ref sequence
            << std::setw(6) << std::right << std::fixed << std::setprecision(6) << swScore  // Smith Waterman score
            << "   " << cigar             // cigar score
//            << " " << ref_ID_Table.getValue(k)                         // Sequence ID (from fasta file
            << std::endl;
        
        

        outFile << oss.str();
    }
    // Close the file
    outFile.close();
    std::cout << "Number of matches : " << resultsSW.size() << std::endl;
    std::cout << "Number of unique read matches : " << unique << std::endl;
    std::cout << "Content written to file: " << newFilename << std::endl;
}



/*
     Check if the rootDir is a regular file
 */
bool isRegularFile(const std::string& filePath) {
    struct stat fileInfo;
    return (stat(filePath.c_str(), &fileInfo) == 0) && S_ISREG(fileInfo.st_mode);
}


/*
     GET FILE NAMES
 */
std::vector<std::string> getFileNames(const std::string& rootDir) {
    std::vector<std::string> fileNames;

    // Check if the rootDir is a regular file and matches the regular expression
    if (isRegularFile(rootDir))
    {
        std::string baseFilename = rootDir;
        std::size_t lastSlashPos = baseFilename.find_last_of("/\\");
        if (lastSlashPos != std::string::npos) {
        baseFilename = baseFilename.substr(lastSlashPos + 1);
        }
        fileNames.push_back(rootDir);
        return fileNames;
    }
    // Open the directory
    DIR* directory = opendir(rootDir.c_str());
    if (directory == nullptr) {
        // Handle error if failed to open directory
        return fileNames;
    }
    // Read directory entries
    dirent* entry;
    while ((entry = readdir(directory)) != nullptr) {
        std::string fileName = entry->d_name;
        std::string filePath = rootDir + "/" + fileName;

        if (fileName == "." || fileName == "..") {
            // Skip current directory and parent directory entries
            continue;
        }

        if (isRegularFile(filePath)) {
            fileNames.push_back(filePath);
        } else {
            // Recursively process subdirectories
            std::vector<std::string> subDirFiles = getFileNames(filePath);
            fileNames.insert(fileNames.end(), subDirFiles.begin(), subDirFiles.end());
        }
    }

    // Close the directory
    closedir(directory);

    return fileNames;
}

/*
    EXTRACT REF NAME
 */
std::string extractRefName(const std::string& filePath) {
    // Extract the file name from the file path
    std::string fileName = filePath.substr(filePath.find_last_of('/') + 1);

    // Remove the file extension
    std::string baseName = fileName.substr(0, fileName.find_last_of('.'));

    return baseName;
}


/*
    GET BASE NAME
 */
std::string getBaseName(const std::string& filePath) {
    std::stringstream ss(filePath);
    std::string baseName;
    std::getline(ss, baseName, '.');
    return baseName;
}



/*
    GET REFERENCE
*/
void getReference(std::string inputFile,
                  std::vector<uint32_t>& innerMapBlob,
                  std::vector<uint32_t>& outerMapBlob,
                  std::string& genomeStr,
                  uint32_t& MASK,
                  QueryTable& ref_ID_Table)
{
    
    std::string outerMapFilename = getBaseName(inputFile)+"_"+std::to_string(KMERSIZE)+"_OuterBlob.idx";
    std::string innerMapFilename = getBaseName(inputFile)+"_"+std::to_string(KMERSIZE)+"_InnerBlob.idx";
    std::string genomeFilename = getBaseName(inputFile)+"_genome.idx";
    std::string refIDFilename = getBaseName(inputFile)+"_refID.idx";

    // load ref ID table
    ref_ID_Table = readTableFromFile(refIDFilename);
    
    // DeSerialize the map
    auto start = std::chrono::steady_clock::now();
    genomeStr = readTextBlobFromFile(genomeFilename);
    deserializeMap(innerMapFilename, outerMapFilename, innerMapBlob, outerMapBlob);
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    int minutes = (int) duration.count() / (1000 * 60);
    float seconds = ((duration.count() - minutes * 1000 * (float) 60))/1000;
    std::cout << "DeSerialising time: " << minutes << " min " << seconds << " sec" << std::endl;

    int numBitsToKeep = std::ceil(std::log2(genomeStr.size()));
    if (numBitsToKeep == 32) {
        MASK = UINT32_MAX; // Set all bits to 1
    } else {
        MASK = (1 << numBitsToKeep) - 1;
    }
//    std::cout << "MASK: " << std::hex << MASK << std::dec << std::endl;
}


/*
     LOAD SAMPLES
 */
void loadSamples(const std::string& fastaFile) {
    // Read sample file into memory
    std::cout << std::endl << "Loading Samples: " << fastaFile << std::endl;
    std::string line;
    std::string sequence; // Variable to collect the DNA sequence
    uint64_t fileSize;

    // Start the timer
    auto start = std::chrono::steady_clock::now();
    char *data = read_entire_file(fastaFile.c_str(),fileSize);
    if (data==NULL) {
        cerr << "Failed to read " << fastaFile << endl;
        exit(1);
    }
    long minLen= (std::numeric_limits<long>::max)();
    long maxLen= 0;

    char *p=data-1; // starrt of line
    char *q; // end of line
    do
    {
         p++;
         q = strchr(p, '\n'); // find end of line
         if (*p=='>') {
             if (!sequence.empty()) {
                 sampleSequences.push_back(sequence);
                 if (sequence.size() < minLen) minLen= (int) sequence.size();
                 if (sequence.size() > maxLen) maxLen= (int) sequence.size();
                 sequence.clear(); // Clear the sequence for the next DNA sequence
             }
             // now skip sequence info
             p=q;
             continue;
         }
         if (q==NULL) {
             // hit last line
             q = strchr(p,'\0'); // find end of data
             line.assign(p,q);
             sequence += line;
             sampleSequences.push_back(sequence); // push last sequence
             break;
         }
         line.assign(p,q);
         sequence += line;
         p=q;
     } while (p != NULL && *(p + 1) != '\0');
     if (!sequence.empty()) {
         sampleSequences.push_back(sequence);
         if (sequence.size() < minLen) minLen= (int) sequence.size();
         if (sequence.size() > maxLen) maxLen= (int) sequence.size();
     }
    free(data);
    // Stop the timer
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    int minutes = (int) duration.count() / (1000*60);
    int seconds = (duration.count()/1000) % 60;
    std::cout << sampleSequences.size() << " samples, loading time: "  << minutes << " min " << seconds << " sec" << std::endl;
    std::cout << "  min seq length " << minLen << ", max seq length " << maxLen << std::endl;
}


/*
    FIND MATCHES
 */
void findMatches(std::string genomeStr, int minMatches, int skip, int startIndex, int endIndex, std::vector<uint32_t>& innerMapBlob, std::vector<uint32_t>& outerMapBlob, uint32_t MASK, int chunkSize) {
    // Start the timer
    // auto start = std::chrono::steady_clock::now();

    // localResults stores the sample index and matching refID, minPos, maxPos of seeds found
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> localResults;
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, double, std::string>> localResultsSW;
    std::set<uint32_t> hitSet; // for keeping positions of seed hits on reference
    uint32_t kmerHash;

    for (int32_t i = startIndex; i <= endIndex; ++i)
    {
    //    int z=95; // debugging
    //    for (uint32_t i = z; i <= z; ++i)
    //    {
    //        if(i==z) {
    //            std::cout << "i = " << i << std::endl;
    //            std::cout << sampleSequences[i] << std::endl;
    //            std::string ref = genomeStr.substr(137105276,151);
    //            std::cout << ref << std::endl;
    //        }

        // process sample read
        hitSet.clear(); //Clear the set of seed hits.
        int sequenceLength = (int) sampleSequences[i].size(); //Get the length of the current sample sequence.
        int lastPos = sequenceLength-KMERSIZE; // last seed position not exceeding read length
        int numSeeds = ceil(lastPos/(double) skip)+1; // number of seeds needed to cover the read
        std::vector<std::set<uint32_t>> inputSets(numSeeds);// A vector of sets to store positions of seed hits for each seed
        int seedMatches=0; // Count of seed matches.
        int maxMisses = numSeeds - MIN_MATCHES; // maximum allowed misses before skipping the rest of the seeds.
        std::string prevSeed; //String variable to store the previous seed
        int startPos; // The starting position for the current seed
 
        //std::cout << "sequence length: " << sequenceLength << std::endl;
        //std::cout << "maxMisses: " << maxMisses << std::endl;


        for (int seedIndex=0; seedIndex < numSeeds; seedIndex++) {
            startPos = seedIndex*skip; // move seed position along
            startPos = (startPos > lastPos) ? lastPos : startPos; // avoid overshoot of last seed

            // Extract a seed (a kmer) from the sample sequence
std::cout << startPos << ":" << KMERSIZE << "\n";
std::cout << sampleSequences[i] << "\n";
            std::string seed = sampleSequences[i].substr(startPos, KMERSIZE);
std::cout << "seed:" << seed << "\n";
            uint64_t pkmer = packKmer(seed.c_str());
            pkmer = pkmer ^ reverse_complement(pkmer, KMERSIZE);// canonical kmer

            if (MASK<UINT32_MAX)
                kmerHash = murmurHash3(pkmer) & MASK;
            else
                kmerHash = murmurHash3(pkmer);

            // Lookup the seed matches in the kMerMap
            std::vector<uint32_t> innerVector = getInnerVector(innerMapBlob, outerMapBlob, kmerHash);
std::cout << "innerVector empty:" << innerVector.empty() << "\n";
             if (!innerVector.empty()) {
                 seedMatches++; // here if seed has matches (positions in the genome)
                 inputSets[seedIndex].insert(innerVector.begin(), innerVector.end()); // insert matches to current seed hits
             }
             else {
                 if ((seedIndex - seedMatches)>maxMisses)
                     break;
             }
        }

std::cout << "Matches:" << seedMatches << "\n";
         if (seedMatches<MIN_MATCHES)
            continue;// read does not meet minimum seed mathces threshold
        
        std::vector<std::pair<uint32_t, uint32_t>> validSets = validSpans(inputSets, minMatches, lastPos);
        
        for (const auto& validSet : validSets) {
            localResults.push_back(std::make_tuple(i, validSet.first, validSet.first, validSet.second));
        }
    }
    
    // std::cout << "localResults: " << localResults.size() << std::endl;
    // bool CALCSW = true; // flag indicating if SW is required

    if (GET_SW) {
        if (localResults.size() >0) {
            // set Smith Waterman scoring matrix.
            // this setting gives the number of base matches for best alignment.
            int matchScore = 1;
            int mismatchScore = -1;
            int gapScore = -1;
            localResultsSW = getSW(genomeStr, localResults, matchScore, mismatchScore, gapScore, CUTOFF);
        }
    }
    else
    {
        // here if no SW required
        for (const auto& result : localResults) {
            uint32_t sampleIndex, referenceIndex, minPos, maxPos;
            std::tie(sampleIndex, referenceIndex, minPos, maxPos) = result;
            std::string noCigar = "*";
            localResultsSW.push_back(std::make_tuple(sampleIndex, minPos, maxPos, CUTOFF, noCigar));
        }
    }
    
    // std::cout << ((double) startIndex)/chunkSize << " After SW: " << localResultsSW.size() << std::endl;
    // Merge local results into global results vector
    if (localResultsSW.size() >0) {
        std::lock_guard<std::mutex> lock(resultsMutex);
        resultsSW.insert(resultsSW.end(), localResultsSW.begin(), localResultsSW.end());
    }
//    std::cout << ((double) startIndex)/chunkSize <<  " maxList: " << maxList << std::endl;
//    auto end = std::chrono::steady_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//    int minutes = (int) duration.count() / (1000 * 60);
//    float seconds = ((duration.count() - minutes * 1000 * (float) 60))/1000;
//    std::cout << endl << ((double) startIndex)/chunkSize <<  " " << "Total time: " << minutes << " min " << seconds << " sec" << std::endl;
}


/*
    PAR FIND MATCHES
 */
void parFindMatches(int N, int M, int numThreads,
                    std::vector<uint32_t>& innerMapBlob,
                    std::vector<uint32_t>& outerMapBlob,
                    std::string genomeStr,
                    uint32_t MASK)
{
    auto start = std::chrono::steady_clock::now();
    std::vector<std::thread> threads;
    int chunkSize = (int) sampleSequences.size() / numThreads;
//    std::cout << "Chunk Size " << chunkSize << std::endl;
    int startIndex = 0;
    int endIndex = 0;
     for (int i = 0; i < numThreads; ++i) {
        startIndex = i * chunkSize;
        endIndex = (i == numThreads - 1) ? (int) sampleSequences.size() - 1 : (int) (startIndex + chunkSize - 1);

         threads.emplace_back(findMatches, genomeStr, N, M, startIndex, endIndex, std::ref(innerMapBlob), std::ref(outerMapBlob), MASK, chunkSize);
    }

    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }
    // Sort the resultsSW vector in descending order based on the SW scores
    std::sort(resultsSW.begin(), resultsSW.end(), SortByScore());

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    int minutes = (int) duration.count() / (1000 * 60);
    float seconds = ((duration.count() - minutes * 1000 * (float) 60))/1000;
//    std::cout << endl << "Found " << resultsSW.size() << " matches" << std::endl;
    std::cout << std::endl << "Matching time: " << minutes << " min " << seconds << " sec" << std::endl;
}


/*
    INITIALISE KISS
 */

void intialiseKISS(int argc, char* argv[]) {
    std::string arg;
    if (argc>1) {
        arg = argv[1];
    }
    if ((argc==1) || arg=="-help") {
        cout << "KISS usage:  ./KISS -argName1 -argValue1 -argName2 -argValue2 ... " << endl;
        cout << "ArgMames are as follows:" << endl;
        cout << "-match_all indicate whether all matches (true|false)" << endl;
        cout << "-sw indicates whether Smith Waterman alignment score is required (true|false)" << endl;
        cout << "-kmerSize is the number of bases in the seed kmer (recommended 16 to 32)" << endl;
        cout << "-min_matches is the minimum number of seed matches in the reference" << endl;
        cout << "-seed_skip is the number of bases skipped for successive seeds" << endl;
        cout << "-cut_off is the minimal SW score to report matches for" << endl;
        cout << "-threads is the number of parallel threads (defeult is HW concurrency)" << endl;
        cout << "-refernce is the name of the reference file" << endl;
        cout << "-data is the name of the data file to be filtered" << endl;
        cout << "-out_dir is the name of directory where results go" << endl;
        cout << "-regex is for recursive search of data files" << endl;
        cout << "example ./filter -kmer_size 30 -min_matches 3 -seed_skip 30 -cut_off 100 -reference CutibacteriumGenome.fasta -data gCSD16_NYC_17_1.fasta -out_dir temp" <<endl;
        exit(0);
    }

    // some default values (overide with cmd line arguments, see below)
    KMERSIZE = 32;
    MIN_MATCHES = 5;
    SEED_SKIP = 256;
    THREADS = std::thread::hardware_concurrency();
    CUTOFF = 100;
    REFERENCE = "";/* "/Users/geva/Crispr/AMR106.fasta"; */
    READS = ""; /* /Users/geva/CAMDA/MatlabCamda2023/gCSD16_NYC_17_1.fasta" */;
    OUTPUT_DIR = "save";
    std::string regexPattern = "";  // Example: Match all .fasta files with ".*\\.fasta";
    MATCH_ALL = false;
    GET_SW = true;

    for (int i = 1; i < argc; ++i) {
        arg = argv[i];
        if (i + 1 >= argc) {
            std::cout << "Error: Missing value for " << arg << " option." << std::endl;
            continue;
        }

        std::string value = argv[i + 1];
        ++i;  // Skip the next argument

        // Process the command line option
        if (arg == "-kmer_size") {
            KMERSIZE = std::stoi(value);
        } else if (arg == "-min_matches") {
            MIN_MATCHES = std::stoi(value);
        } else if (arg == "-seed_skip") {
            SEED_SKIP = std::stoi(value);
        } else if (arg == "-cut_off") {
            CUTOFF = std::stoi(value);
        } else if (arg == "-threads") {
            THREADS = std::stoi(value);
        } else if (arg == "-reference") {
            REFERENCE = value;
        } else if (arg == "-data") {
            READS = value;
        } else if (arg == "-out_dir") {
            OUTPUT_DIR = value;
        } else if (arg == "-match_all") {
            MATCH_ALL = (value=="true");
        } else if (arg == "-sw") {
            GET_SW = (value=="true");
        } else {
            cerr << "Error: Unknown option: " << arg << std::endl;
        }
    }
    // set up gobal encodingTable for packKmer
    kmer_encoding_table[(size_t)'A'] = 0;
    kmer_encoding_table[(size_t)'C'] = 1;
    kmer_encoding_table[(size_t)'G'] = 2;
    kmer_encoding_table[(size_t)'T'] = 3;
    
    std::string rootDir = READS;
    fileNames = getFileNames(rootDir);
    // Print the file names
    std::cout << "Processing " << fileNames.size() << " files" << std::endl;
        //    for (const auto& fileName : fileNames) {
        //        std::cout << fileName << std::endl;
        //    }
    
    //ESZTERQ: why be able to set min_matches if the change is over run anyway? 
    MIN_MATCHES = ceil((CUTOFF - KMERSIZE)/(double) SEED_SKIP)+1;

    cout << "KISS run parameters" << endl;
    cout << "kmer_size: " << KMERSIZE << endl;
    cout << "min_matches: " << MIN_MATCHES << endl;
    cout << "seed_skip: " << SEED_SKIP << endl;
    cout << "threads: " << THREADS << endl;
    cout << "cut_off: " << CUTOFF << endl;
    cout << "reference: " << REFERENCE << endl;
    cout << "data: " << READS << endl;
    cout << "output directory: " << OUTPUT_DIR << endl << endl;
    cout << "      find all matches: " << (MATCH_ALL==true ? "true": "false") << endl;
    cout << "compute Smith Waterman: " << (GET_SW==true ? "true": "false") << endl << endl;

    if (REFERENCE.empty() || READS.empty()) {
        cerr << "Reference file name and Data file name must be provided" << endl;
        exit(1);
    }
}
/*
     MAIN
 */
int main(int argc, char *argv[]) {

    std::vector<uint32_t> innerMapBlob;
    std::vector<uint32_t> outerMapBlob;
    std::string genomeStr;
    uint32_t MASK;
    QueryTable ref_ID_Table;
    
    // Start the timer
    auto startAll = std::chrono::steady_clock::now();

    // set up KISS parameters
    intialiseKISS(argc, argv);

    getReference(REFERENCE,innerMapBlob,outerMapBlob,genomeStr,MASK,ref_ID_Table); // load the reference index

    size_t startIndex = 0; // Specify the index from which to start (in case restarting after crash)
    for (size_t i = startIndex; i < fileNames.size(); ++i) {
        const std::string& filename = fileNames[i];
        std::cout << i << " ******** processing " << filename << " ********" << std::endl;
        // Process the current file
        // Start the timer
        resultsSW.clear();
        sampleSequences.clear();
        loadSamples(filename); // load the reads (to be filtered against the reference)
        parFindMatches(MIN_MATCHES, SEED_SKIP, THREADS,innerMapBlob,outerMapBlob,genomeStr,MASK); // find read matches in the refrence
        writeMatchesToFile(filename,resultsSW,ref_ID_Table);
    }
    // report overall program duration
    auto endAll = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endAll - startAll);
    int minutes = (int) duration.count() / (1000 * 60);
    float seconds = ((duration.count() - minutes * 1000 * (float) 60))/1000;
    std::cout << endl << "Total time: " << minutes << " min " << seconds << " sec" << std::endl;

    return 0;
}

