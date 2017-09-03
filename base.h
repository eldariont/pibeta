#ifndef SEQAN_HEADER_PACMAPPER_H
#define SEQAN_HEADER_PACMAPPER_H

using namespace seqan;

//===================================================================
// const var and type def
//===================================================================
//template <typename TSpec = void>

struct Const_{
    
    typedef unsigned _INT_;
    typedef uint64_t _LLT_;
    typedef float _FLT_;
    typedef CharString _CSR_;
    typedef _CSR_   _FILE_;
    typedef _INT_   _SHAPE_;
    typedef seqan::Dna5 _CA5_;

    static const _SHAPE_ _SHAPELEN;
    static const _SHAPE_ _SHAPEWHT;
    static const _INT_ _BLOCKSIZE;
    static const _INT_ _DELAT; 
    static const _INT_ _THRESHOLD; 
    static const _FLT_ _ALPHA ;
    static const _INT_ _KMERSTEP;
        
//    void operator (_LLT_ & a, _LLT_ & b)
};

const typename Const_::_FLT_ Const_::_ALPHA = 0.8;
const typename Const_::_SHAPE_ Const_::_SHAPELEN = 25;
const typename Const_::_SHAPE_ Const_::_SHAPEWHT = 18;
const typename Const_::_INT_ Const_::_BLOCKSIZE =  1000;
const typename Const_::_INT_ Const_::_DELAT = 32; 
const typename Const_::_INT_ Const_::_THRESHOLD = 30; 
const typename Const_::_INT_ Const_::_KMERSTEP = 1000;
     
struct Options{
    typename Const_::_INT_  kmerLen;
    typename Const_::_INT_  MiKmLen;
    typename Const_::_FILE_ rPath;
    typename Const_::_FILE_ gPath;
    typename Const_::_FILE_ oPath;
    bool        Sensitive; 
    Options():
        kmerLen(Const_::_SHAPELEN),
        MiKmLen(Const_::_SHAPEWHT),
        rPath(""),
        gPath(""),
        oPath(""),
        Sensitive(false)
        {}
} options;
//
//template <typename TDna>
//struct RecordBase
//{
//    typedef typename Const_::_CSR_ RecId;
//    typedef String<TDna> RecSeq; 
//};
//
//template<typename TDna>
//struct PMRecord
//{
//    typedef RecordBase::RecId RecId;
//    typedef RecordBase::RecSeq RecSeq;
//    typedef String<RecId> RecId;
//    typedef StringSet<RecSeq> RecSeq;
//
//    RecId id;
//    RecSeq seq1, seq2; //seq1=read, seq2=ref
//    int loadRecord(typename Const_::_FILE_ const & path);
//};
//
//struct AnchorBase{
//    typedef typename Const_::_INT_ AnchorType_; 
//    AnchorType_ anchor;    
//
//    AnchorType_ setAnchorNode(AnchorType_ anchorPos,  AnchorType_ kmerPos);
//    AnchorType_ getAnchorPos();
//    AnchorType_ getKmerPos();
//};
//
//struct Anchor{
//    static const Const_::_INT_ size = 131072;
//    SimpleAnchor_ set[size];
//    Anchor(){
//        //set = {1ULL << 63};
//    };
//    void init();
//    SimpleAnchor_ & operator [](unsigned);
//};
//
//
//template <typename TDna = typename Const_::_CA5_, typename CoreMinimizer>
//struct CoreBase{
//    typedef typename Const_::_CA5_ CoreDefault
//    typedef typename Const_::_INT_ ShapePar;
//    typedef typename PMReocrd<TDna>::RecSeq RecSeq; 
//    typedef Shape<TDna, CoreMinimizer> CoreShape;
//    typedef Index<RecSeq, IndexQGram<CoreMinimizer, OpenAddressing> > CoreIndex;
//
//};
//
//template <typename TDna = CoreBase<>::CoreDefault, typename Minimizer>
//struct PMCore
//{
//    typedef typename CoreBase<TDna, Minimizer>::CoreIndex Index;
//    typedef typename CoreBase<TDna, Minimizer>::RecSeq Seq;
//    Index index;
//    Anchor anchor;
//}
//
//template <unsigned bit=32,uint64_t n1, uint64_t n2>
//struct makeInt{
//    enum{Value = (n1<< bit) + n2}
//}
//
//template <unsigned bit=32,uint64_t n1, uint64_t n2>
//struct getInt{
//    enum{Value = (n1<< bit) + n2}
//}
//
////template <typename TSpec = void>
//struct ResBase{
//    typedef typename Const_::_LLT_ SeqId;
//    typedef typename Const_::_INT_ SeqLen;
//    typedef typename Const_::_LLT_ MapPos;
//    typedef typename Const_::_LLT_ MapStrand;
//    typedef typename Const_::_LLT_ MapScore;
//
//    static const Const::_INT_ bit = 32
//    static const Const::_LLT_ mask = (1 << idBit) - 1;
//    MapStrand   strand = 1ULL << 63;
//
//};
//
////template<typename TSpec = void>
//struct PMRes{
//    typedef typename ResBase::SeqId Id;
//    typedef typename ResBase::MapPos Pos;
//    typedef typename ResBase::MapScore Score;
//    typedef typename ResBase::MapStrand Strand;
//    typedef String<Id> ResID;
//    typedef String<Pos> ResPos;
//    typedef String<Score> ResScore;
//    
//    ResID   id;
//    ResPos  p1, p2; //p1=start, p2=end
//    ResScore score;  
//    MapStrand strand;
//
//    SeqId getId1(SeqId id); 
//    SeqId getId2(SeqId id);
//    MapPos getPos1(MapPos pos);
//    MapPos getPos2(MapPos pos);
//    createId(Id id1, Id id2)
//    createPos(Pos p1, Pos p2, Strand strand)
//
//    MapStrand getStrand();
//    void operator ()(SeqId & id, MapStrand & strand, MapPos& pos);
//};
//
//inline PMRes::SeqId PMRes::getId1(PMRes::SeqId id)
//{
//    return id >> ResBase::bit;
//}
//
//inline PMRes::MapPos PMRes::getPos2(PMRes::MapPos p)
//{
//    return p & ResBase::mask;
//}
//inline PMRes::MapPos PMRes::getPos1(PMRes::MapPos p)
//{
//    return p >> ResBase::bit;
//}
//
//inline PMRes::MapPos PMRes::getPos2(PMRes::MapPos p)
//{
//    return p & ResBase::mask;
//}
//
//template <typename TDna = typename Const_::_CA5_, typename TSpec = void>
//struct MapperBase
//{
//    PMRecord record;
//    PMCore  core;
//    PMRes   res;
//    void loadRecord(path)
//     
//};
//
//
////template <typename TSpec = void>
struct MapParm{
    typename Const_::_INT_  blockSize;
    typename Const_::_FLT_  alpha;  
    typename Const_::_INT_  delta;
    typename Const_::_INT_  threshold;
    typename Const_::_INT_  kmerStep;
    
    MapParm():
        blockSize(Const_::_BLOCKSIZE),
        alpha(Const_::_ALPHA),
        delta(Const_::_DELAT),
        threshold(Const_::_THRESHOLD),
        kmerStep(Const_::_KMERSTEP)
        {}
    void setMapParm(Options & options);
} _DefaultMapParm;
//
//
//template <typename TDna = typename Const_::_CA5_, typename TSpec = void>
//struct Mapper
//{
//    typename MapperBase<TDna>::MSeqSet reads;
//    typename MapperBase<TDna>::MSeqSet gnome;
//    typename MapperBase<TDna>::MShape shape;
//    typename MapperBase<TDna>::MIndex index;
//
//    MapParm        mapParm;
//    typename MapperStruct<TDna>::MRes res;
//    Mapper()
//    {};
//    //Mapper(Options const & options)
//    //{
//    //    loadRecords(options);
//    //    setMapParm(options);
//    //};
//};
//
////===================================================================
//// function  
////===================================================================
//
//inline typename PMRes::SeqId PMRes::getRdId(unsigned k){
//    //return this->id[k] >> this->_bit;
//    return 0;
//} 
//
//inline typename PMRes::SeqId PMRes::getRfId(unsigned k){
// //   return this->id[k] & this->_mask;
//    return 0;
//}
//
//inline void PMRes::operator() (typename ResRecord_::SeqId_ & mid,
//                typename ResRecord_::MapStrand_ & mstrand,
//                //typename ResRecord_::SeqLen & mlen,
//                typename ResRecord_::MapPos_ & mpos)
//{
//    //appendValue(id, mid);
//    //appendValue(pos, mpos | mstrand);
//    //return 0;
//}
//
//int _output()
//{
//    
//    return 0;
//}
//
///*
//template <typename Stream, typename Records>
//inline void _writeRecords(Stream & stream, Records & records){
//      
//}
//
//tempalte <typename WOptions_, typename Res>
//inline void writeRes(WOptoins & wopt, Res & res)
//{
//    
//}
//
//*/
//
//
//
//template <typename TDna>
//inline int SimpleRecord_<TDna>::loadRecord(typename Const_::_FILE_ const & path)
//{
//    //std::cerr << "loading sequence " << std::endl;
//    seqan::readRecords(this.ids, this.seq, path);
//    return 0; 
//}
///*
//void loadSequence(Options & options, id)
//{
//    try{
//        if (!_loadSeq())
//            throw  
//    }
//}
//*/
//
///*
//struct _compltStr
//{
//     
//    static String<Dna5> _complt = "tgcan";
//    void operator()
//    {
//        resize(res, length(str));
//    for (unsigned k = 0; k < length(str); k++)
//     //   res[k]=_complt[str[k] - 'A'];
//        res[k] = _complt[(unsigned)ordValue(str[k])];
//
//    }
//}
//*/



const unsigned shapelength = 25; 
const unsigned shapeweight = 18; 
//const unsigned blocklimit = 32;

typedef Iterator<String<Dna> >::Type TIter;
typedef Iterator<String<Dna5> >::Type TIter5;
typedef Shape<Dna, Minimizer<shapelength> > TShape;
typedef Shape<Dna5, Minimizer<shapelength> > TShape5;
typedef Shape<Dna, UngappedShape<shapelength> > TShape_u;
typedef Shape<Dna, SimpleMShape> TMShape;

typedef Index<StringSet<DnaString>, IndexQGram<Minimizer<shapelength>, OpenAddressing > > TIndex;
typedef Index<StringSet<String<Dna5> >, IndexQGram<Minimizer<shapelength>, OpenAddressing > > TIndex5;
typedef Index<StringSet<DnaString>, IndexQGram<UngappedShape<shapelength>, OpenAddressing > > TIndex_u;

static String<Dna5> _complt = "tgcan";
inline void _compltStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
     //   res[k]=_complt[str[k] - 'A'];
        res[k] = _complt[(unsigned)ordValue(str[k])];
}

inline void _compltRvseStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
     //   res[k]=_complt[str[k] - 'A'];
        res[k] = _complt[(unsigned)ordValue(str[length(str) - k])];
}


//encapsulate
/*
inline void Anchor::init(){
   // for (unsigned k = 0; k < this->size; k++)
   // {
   //     this->set[k] = 0;
   // } 
}

inline SimpleAnchor_ & Anchor::operator[](unsigned k)
{
    return set[k];
}

inline void typename SimpleAnchor_::setAnchorNode
(typename SimpleAnchor_::AnchorType_ & anchorPos, typename SimpleAnchor_::AnchorType_ & kmerPos)
{
    anchor = (anchorPos << Const_::AnchorPos) + kmerPos;
}

inline typename SimpleAnchor_::AnchorType_ SimpleAnchor_::getAnchorPos();
{
    return anchor >> anchor_bit; 
}

inline typename SimpleAnchor_::AnchorType_ SimpleAnchor_::getKmerPos();
{
    return anchor & anchor_mask;
}

*/

/*
template <TDna>
inline unsigned _mnMapRead(typename Mapper<TDna>::MIndex & index, 
typename Mapper<TDna>::MSeq & read, Anchor & anchor, typename Mapper::MParm & mapParm )
{
    uint64_t x = 1;
    unsigned _block = (mapParm.blockSize < length(read))?mapParm.blockSize:length(read);
    unsigned _dt = _block * (mapParm.alpha / (1 - mapParm.alpha));

    anchor.init();
    for (unsigned h=0; h <= length(read) - _block; h += _dt)
    {
        hashInit(index.shape, begin(read) + h);
        for (unsigned k = h; k < h + _block; k++)
        {
            hashNext(index.shape, begin(read) + k);
            uint64_t dn = getDir(index, index.shape);
            uint64_t pre = ~0;
            if(_getBodyCounth(index.dir[dn+1]) - _getBodyCounth(index.dir[dn]) < mapParm.delta)
            {
                for (uint64_t n = _getBodyCounth(index.dir[dn]); n < _getBodyCounth(index.dir[dn + 1]); n++)
                {
                    if (index.sa[n] - pre > mapParm.kmerStep)
                    {
                        anchor[x++].setAnchorNode(((index.sa[n]- k) << 20) + k);
                        pre = index.sa[n];
                    }
                }
            }
        }
    }
    anchor[0] = anchor[1];
    uint64_t max = 0, c_b=0, ak=anchor[0], cbb=0, mask_cb = (1<<20) - 1, sb=0;
    unsigned start=0;

    std::sort(begin(anchor), begin(anchor) + x);
    
    for (uint64_t k = 1; k <= x; k++)
    {
        if ()
        if (anchor[k]-ak < (1000<<20))
        {
            cbb++;
                    }
        else
        {
            std::sort(begin(anchor)+sb, begin(anchor)+k, 
            [&mask_cb](uint64_t &a, uint64_t &b){return (a & mask_cb) < (b & mask_cb);});
            for (uint64_t m = sb+1; m < k; m++)
            {
                if(((anchor[m]-anchor[m-1]) & mask_cb) > shapelength) 
                    c_b += shapelength;
                else
                {
                    c_b += (anchor[m] - anchor[m-1]) & mask_cb; 
                }
            }
            if (c_b > max)
            {
                max = c_b;
                start = sb;
            }
            sb = k;
            ak = anchor[k];
            cbb = 1;
            c_b = shapelength;
        }

    }
    return start + (max << 20) ;
}
//encapsulate end
//


*/


template <typename TIndex, typename TObj>
inline unsigned _mnMapReads(TIndex & index, String<TObj> & read,  uint64_t* const anchor, MapParm &  mapParm)
{
    hashInit(index.shape, begin(read));
    for (uint64_t h = 0; h < length(anchor); h++)
        anchor[h] = 0;
    //$ anchor.init();
    uint64_t x = 1;
    //uint64_t block =(1000 < length(read))?1000:length(read);
    //uint64_t l = block * (0.8/0.2);

    unsigned _block = (mapParm.blockSize < length(read))?mapParm.blockSize:length(read);
    unsigned _dt = _block * (mapParm.alpha / (1 - mapParm.alpha));

    for (unsigned h=0; h <= length(read) - _block; h += _dt)
    //for (unsigned h=0; h <= length(read) - _block; h += _block)
    {
        hashInit(index.shape, begin(read) + h);
        for (unsigned k = h; k < h + _block; k++)
        {
            hashNext(index.shape, begin(read) + k);
            uint64_t dn = getDir(index, index.shape);
            uint64_t pre = ~0;
            //if(_getBodyCounth(index.dir[dn+1]) - _getBodyCounth(index.dir[dn]) < 32)
            if(_getBodyCounth(index.dir[dn+1]) - _getBodyCounth(index.dir[dn]) < mapParm.delta)
            {
                for (uint64_t n = _getBodyCounth(index.dir[dn]); n < _getBodyCounth(index.dir[dn + 1]); n++)
                {
                    //if (index.sa[n] - pre > 1000)
                    if (index.sa[n] - pre > mapParm.kmerStep)
                    {
                        anchor[x++] = (((index.sa[n])- k) << 20) + k;
                        //anchor.anchor[x++] = anchor.makeAnchorNode(index.sa[n] - k,  k);
                        pre = index.sa[n];
                    }
                }
            }
        }
    }
    anchor[0] = anchor[1];
    uint64_t max = 0, c_b=0, ak=anchor[0], cbb=0, mask_cb = (1<<20) - 1, sb=0;
    unsigned start=0;

    std::sort(begin(anchor), begin(anchor) + x);
    
    for (uint64_t k = 1; k <= x; k++)
    {
        if (anchor[k]-ak < (1000<<20))
        {
            cbb++;
                    }
        else
        {
            std::sort(begin(anchor)+sb, begin(anchor)+k, 
            [&mask_cb](uint64_t &a, uint64_t &b){return (a & mask_cb) < (b & mask_cb);});
            for (uint64_t m = sb+1; m < k; m++)
            {
                if(((anchor[m]-anchor[m-1]) & mask_cb) > shapelength) 
                    c_b += shapelength;
                else
                {
                    c_b += (anchor[m] - anchor[m-1]) & mask_cb; 
                }
            }
            if (c_b > max)
            {
                max = c_b;
                start = sb;
            }
            sb = k;
            ak = anchor[k];
            cbb = 1;
            c_b = shapelength;
        }

    }
    //return anchor[start] >> 20;
    return start + (max << 20) ;
}
/*
template <TDna>
inline void mnMap(Mapper<TDna> & mapper)
{
    typename Mapper<TDna>::MShape shape;
    typename Mapper<TDna>::MIndex index(mapper.genome);
    
    double time = sysTime();
    unsigned res = 0, tmp = 0; 
    Anchor anchor;
    SimpleAnchor_ res;
    
    std::cerr << "Creating index \n";
    createQGramIndexDirOnly(index);
    std::cerr << "Filtreing reads \n";
    
    for (unsigned j = 0; j < length(mapper.reads); j++)
    {   
        res = _mnMapRead(index, mapper.reads[j], anchor, mapper.mapParm);
        if (res.getLength() < mapper.getThreshold())
        {
            _compltRvseStr(mapper.reads[j], rcStr);
            rcRes = _mnMapRead(index, rcStr, anchor, mapper.mapParm);
            res = (rcRes.getLength() > (mapper.getThreshold()) ? rc_res : SimpleAnchor::NULL;           
        }
        //overload ()
        mapper.res(res);
    }
    
}

*/
template <typename TObj>
void mnMap(StringSet<String<TObj> > & genome, StringSet<String<TObj> > & reads, MapParm & mapParm, String<Pair<uint64_t> > & rs)
{
    //typedef ModifiedString<String<Dna5>, ModComplementDna5> Comp;
    TShape5 shape;
    TIndex5 index(genome);
    String<TObj> comStr;
    double time=sysTime();
    unsigned mask1 = (1<<20)-1, res = 0, tmp = 0;
    uint64_t mask  = 131072 - 1;
    uint64_t anchor[mask + 1] = {1ULL<<63};
    String<uint64_t> result;
    resize(result, length(reads));
    std::cerr << "Creating index \n";
    _createQGramIndex(index);
    std::cerr << "Filtering reads\n"; 
 //   _compltStr(reads[0], comStr); 
//    std::cout << reads[0] << std::endl;
//    std::cout << comStr << std::endl;
    for (unsigned j = 0; j< length(reads); j++)
    {
        res = _mnMapReads(index, reads[j], anchor, mapParm);

            //std::cout << (res>>20) << std::endl;
        if (res < (mapParm.threshold << 20))
        {
            _compltRvseStr(reads[j], comStr);
            tmp = _mnMapReads(index, comStr, anchor, mapParm);
            res = (tmp > (mapParm.threshold << 20))?tmp:0;
        }
        
        //appendResult();
        
        //std::cout << j << " " << length(reads[j]) << " " << _getSA_i2(anchor[res & mask1] >> 20) << std::endl;
        //appendValue(rs, makePairanchor[res & mask1] >> 20);
        assignValueI1(rs[j], _getSA_i2(anchor[res & mask1] >> 20));
        assignValueI2(rs[j], _getSA_i2(anchor[res & mask1] >> 20) + length(reads[j]));
        
    }
    std::cerr << sysTime() - time << std::endl;
}

//template <typename TSpec = void>
//void printInfo()
//{
//    
//}

/*
void map(Mapper & mapper)
{
    mnMap(mapper);
}
*/

template <typename TObj>
String<Pair<uint64_t, uint64_t> > map(StringSet<String<TObj> > & genome, StringSet<String<TObj> > & reads, MapParm & mapParm = _DefaultMapParm)
{
    String<Pair<uint64_t, uint64_t> > result;
    resize(result, length(reads));
    //$map.init(genome, reads, mapParm);
    mnMap(genome, reads, mapParm, result);
    //$mnMap(map.ref(), map.reads(), map.mapParm(), map.res();
    for (uint64_t k = 0; k < length(result); k++)
        std::cout << k << " " << result[k].i1 << " " << result[k].i2 << std::endl;
    //$map.printRes();
    return result;
}

#endif
