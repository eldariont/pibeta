

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

    bandGenerator._colOffset = (lastBlock.second - firstBlock.second + 1);
    resize(bandGenerator._gridBlocks,
           (lastBlock.first - firstBlock.first + 1) * bandGenerator._colOffset,
           false,
           Exact{});
    for (auto & block : grid)
    {
        // check that current block is not outside of matrix
        if (firstBlock.second <= block.second <= lastBlock.second && firstBlock.first <= block.first <= lastBlock.first)
        {
            bandGenerator._gridBlocks[block.first * bandGenerator._colOffset + block.second] = true;
        }        
    }

    //TODO(rrahn): Second phase: Identify and connect lose ends with parallelogram
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
