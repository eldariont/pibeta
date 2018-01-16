

#include "base.h"

#include <limits>


struct DPStaticBandGenerator
{
    enum class DIM : bool
    {
        HORIZONTAL,
        VERTICAL
    };
    constexpr static bool HORIZONTAL_DIM{true};
    constexpr static bool VERTICAL_DIM{false};
    // The last block included in the band in horizontal and vertical direction.
    uint32_t horizontalGridIndex{std::numeric_limits<uint32_t>::max()};
    uint32_t verticalGridIndex{std::numeric_limits<uint32_t>::max()};

    String<bool, Packed<>> _gridBlocks{};  // blocks to create the task graph on.
    uint32_t               _colOffset{0};  // offset for the column size

    inline bool _atEnd(size_t const pos, DPStaticBandGenerator::DIM const dim)
    {
        SEQAN_ASSERT(pos < length(_gridBlocks));

        return (dim == DIM::HORIZONTAL) ? ((pos / _colOffset) == ((length(_gridBlocks) - 1) / _colOffset))
                                        : ((pos % _colOffset) == ((length(_gridBlocks) - 1) % _colOffset));
    }

    inline bool _atBegin(size_t const pos, DPStaticBandGenerator::DIM const dim)
    {
        SEQAN_ASSERT(pos < length(_gridBlocks));

        return (dim == DIM::HORIZONTAL) ? ((pos / _colOffset) == 0)
                                        : ((pos % _colOffset) == 0);
    }

    // We need to offer some functionality to access the corresponding coordinates.
    inline bool hasSuccessor(size_t const pos, DPStaticBandGenerator::DIM const dim)
    {
        SEQAN_ASSERT(pos < length(_gridBlocks));

        if (empty(_gridBlocks) || _atEnd(pos, dim))
            return false;

        return (dim == DIM::HORIZONTAL) ? (_gridBlocks[pos + _colOffset]) : (_gridBlocks[pos + 1]);
    }

    inline bool hasPredecessor(size_t const pos, DPStaticBandGenerator::DIM const dim)
    {
        SEQAN_ASSERT(pos < length(_gridBlocks));

        if (empty(_gridBlocks) || _atBegin(pos, dim))
            return false;
        return (dim == DIM::HORIZONTAL) ? (_gridBlocks[pos - _colOffset]) : (_gridBlocks[pos - 1]);
    }
};

template <typename TSeqH, typename TSeqV>
void computeGrid(DPStaticBandGenerator & bandGenerator,
                 TSeqH const & seqH,
                 TSeqV const & seqV,
                 uint32_t const grid_size,
                 bool const strand)
{
    //TODO(rrahn): Write me!
    SEQAN_ASSERT_FAIL("Implement me!");
}

struct PacBioBandGenerator : public DPStaticBandGenerator
{};


void connectInGrid(PacBioBandGenerator & bandGenerator, std::pair<uint32_t, uint32_t> lowerRightBlock, std::pair<uint32_t, uint32_t> upperLeftBlock)
{
    std::cerr << "Connect " << upperLeftBlock.first << "," << upperLeftBlock.second << " with " << lowerRightBlock.first << "," << lowerRightBlock.second << std::endl;
    for (int row = upperLeftBlock.first; row <= lowerRightBlock.first; row++)
    {
        for (int col = upperLeftBlock.second; col <= lowerRightBlock.second; col++)
        {
            bandGenerator._gridBlocks[row * bandGenerator._colOffset + col] = true;
        }
    }
}


bool checkProperGrid(PacBioBandGenerator & bandGenerator)
{
    int cols = bandGenerator._colOffset;
    int rows = length(bandGenerator._gridBlocks) / cols;
    std::cerr << "Grid has " << rows << " rows and " << cols << " cols." << std::endl;

    bool unconnectedBlocksFound = false;

    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            if (bandGenerator._gridBlocks[row * bandGenerator._colOffset + col] == true)
            {
                // first, find out whether current block is connected to block above or left
                bool connected = true;
                if (row == 0 and col == 0)
                {
                    ;
                }
                else if (col == 0)
                {
                    if (bandGenerator._gridBlocks[(row-1) * bandGenerator._colOffset + col] == false)
                    {
                        connected = false;
                    }
                }
                else if (row == 0)
                {
                    if (bandGenerator._gridBlocks[row * bandGenerator._colOffset + (col - 1)] == false)
                    {
                        connected = false;
                    }
                }
                // else
                else{
                    if (bandGenerator._gridBlocks[(row-1) * bandGenerator._colOffset + col] == false && bandGenerator._gridBlocks[row * bandGenerator._colOffset + (col - 1)] == false)
                    {
                        connected = false;
                    }
                }

                if (!connected)
                {
                    std::cerr << "Block " << row << "," << col << " is not connected. Try to connect.." << std::endl;
                    unconnectedBlocksFound = true;
                }

                // now, try to connect unconnected block
                int i = -1;
                while (!connected && i >= std::min(-col, -row))
                {
                    int r = row + i;
                    if (r >= 0)
                    {
                        for (int c = col; c > col + i && c >= 0; c--)
                        {
                            if (bandGenerator._gridBlocks[(r) * bandGenerator._colOffset + c])
                            {
                                connectInGrid(bandGenerator, std::make_pair(row, col), std::make_pair(r, c));
                                connected = true;
                                break;
                            }
                        }
                    }                        
                    if (!connected)
                    {
                        int r = row + i;
                        int c = col + i;
                        if (r >= 0 && c >= 0)
                        {
                            if (bandGenerator._gridBlocks[(r) * bandGenerator._colOffset + c])
                            {
                                connectInGrid(bandGenerator, std::make_pair(row, col), std::make_pair(r, c));
                                connected = true;
                                break;
                            }
                        }
                    }
                    if (!connected)
                    {
                        int c = col + i;
                        if (c >= 0)
                        {
                            for (int r = row; r > row + i && r >= 0; r--)
                            {
                                if (bandGenerator._gridBlocks[(r) * bandGenerator._colOffset + c])
                                {
                                    connectInGrid(bandGenerator, std::make_pair(row, col), std::make_pair(r, c));
                                    connected = true;
                                    break;
                                }
                            }
                        }
                    }
                    i -= 1;
                }

                // if no connection could be made, connect with upper left corner of matrix
                if (!connected)
                {
                    connectInGrid(bandGenerator, std::make_pair(row, col), std::make_pair(0, 0));
                }
            }
        }
    }
    return !unconnectedBlocksFound;
}


//TODO(rrahn): Model band generation.
template <typename TCoordinates>
inline void
computeGrid(PacBioBandGenerator & bandGenerator,
            TCoordinates const & readMappingCoordinates,
            uint32_t const grid_size,
            bool const strand)
{
    std::cerr << "Strand: " << std::boolalpha << strand << std::endl;
    std::vector<std::pair<uint32_t, uint32_t>> grid;
    
    // find lowest position in read and genome
    unsigned genome_start = -1;
    unsigned read_start = -1;
    for (auto && tup : readMappingCoordinates)
    {
        if (std::get<2>(tup) < genome_start)
        {
            genome_start = std::get<2>(tup);
        }
        if (std::get<0>(tup) < read_start)
        {
            read_start = std::get<0>(tup);
        }
    }

    // compute grid indices for each block
    for (auto && tup : readMappingCoordinates)
    {
        auto x_pos = std::get<2>(tup) - genome_start;
        auto y_pos = std::get<0>(tup) - read_start;
        grid.push_back(std::make_pair(x_pos / grid_size, y_pos / grid_size));
        grid.push_back(std::make_pair((x_pos + 192) / grid_size, y_pos / grid_size));
        grid.push_back(std::make_pair(x_pos / grid_size, (y_pos + 192) / grid_size));
        grid.push_back(std::make_pair((x_pos + 192) / grid_size, (y_pos + 192) / grid_size));
        // std::cerr << read << "\t" << mapping << "\t" << std::get<0>(tup) << "\t" << std::get<1>(tup) << "\t" << std::get<2>(tup) << std::endl;
    }
    std::stable_sort(grid.begin(), grid.end());
    grid.erase(std::unique(grid.begin(), grid.end()), grid.end());

    auto & lastBlock  = back(grid);
    auto & firstBlock = front(grid);

    // initialize matrix of blocks to process
    bandGenerator._colOffset = (lastBlock.second - firstBlock.second + 1);
    resize(bandGenerator._gridBlocks,
        (lastBlock.first - firstBlock.first + 1) * bandGenerator._colOffset,
        false,
        Exact{});
    
    // fill in blocks that have to be processed
    for (auto & block : grid)
    {
        // check that current block is not outside of matrix
        if (firstBlock.second <= block.second <= lastBlock.second && firstBlock.first <= block.first <= lastBlock.first)
        {
            bandGenerator._gridBlocks[block.first * bandGenerator._colOffset + block.second] = true;
        }        
    }

    // make sure that all blocks are connected
    checkProperGrid(bandGenerator);
    SEQAN_ASSERT(checkProperGrid(bandGenerator));
}

struct DPVerifierTraits
{
    // The algorithm to choose.
    using TAlgorithmType    = GlobalAlignment_<FreeEndGaps_<True, False, True, False>>;
    // The Gaps to choos
    using TGapType          = AffineGaps;
    // The Band to choose.
    using TBandType         = PacBioBandGenerator;
    // The traceback.
    using TTracebackType    = TracebackOff;
    // The output to choose.
    using TFormat           = ArrayGaps;
};


// template <typename TSpec, typename TSimdSpec,
//           typename TSetH,
//           typename TSetV,
//           typename TSettings,
//           typename TCallable>
// inline void
// alignExecBatch(ExecutionPolicy<WavefrontAlignment<TSpec>, TSimdSpec> const & execPolicy,
//                TSetH const & setH,
//                TSetV const & setV,
//                TSettings const & settings,
//                TCallable && callback)
// {
//     using TSeqH = typename Value<TSetH const>::Type;
//     using TSeqV = typename Value<TSetV const>::Type;
//
// #ifdef SEQAN_SIMD_ENABLED
//     using TExecutor = std::conditional_t<std::is_same<TSimdSpec, Vectorial>::value,
//                                          AsyncWaveAlignExecutorSimd<TSeqH, TSeqV, TSettings, TSpec>,
//                                          AsyncWaveAlignExecutor<TSeqH, TSeqV, TSettings>>;
// #else
//     using TExecutor = AsyncWaveAlignExecutor<TSeqH, TSeqV, TSettings>;
// #endif
//     TExecutor executor(settings, execPolicy);
//
//     for (size_t i = 0u; i < length(setH); ++i)
//     {
//         submit(executor, setH[i], setV[i], std::forward<TCallable>(callback));
//     }
//     wait(executor);
// }
