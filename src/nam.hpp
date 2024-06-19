#ifndef STROBEALIGN_NAM_HPP
#define STROBEALIGN_NAM_HPP

#include <vector>
#include <array>
#include "index.hpp"
#include "randstrobes.hpp"

struct PartialSeed {
    size_t hash;
    unsigned int start;
    bool is_reverse;
    bool operator==(const PartialSeed& rhs) const {
        return (hash == rhs.hash) && (start == rhs.start) && (is_reverse == rhs.is_reverse);
    }
};

struct MainHit {
    int query_start;
    int query_end;
    int ref_start;
    int ref_end;
};

struct RescueHit {
    size_t position;
    unsigned int count;
    unsigned int query_start;
    unsigned int query_end;

    bool operator< (const RescueHit& rhs) const {
        return std::tie(count, query_start, query_end)
               < std::tie(rhs.count, rhs.query_start, rhs.query_end);
    }
};

// Non-overlapping approximate match
struct Nam {
    int nam_id;
    int query_start;
    int query_end;
    int query_prev_hit_startpos;
    int ref_start;
    int ref_end;
    int ref_prev_hit_startpos;
    int n_hits = 0;
    int ref_id;
    float score;
//    unsigned int previous_query_start;
//    unsigned int previous_ref_start;
    bool is_rc = false;

    int ref_span() const {
        return ref_end - ref_start;
    }

    int query_span() const {
        return query_end - query_start;
    }

    int projected_ref_start() const {
        return std::max(0, ref_start - query_start);
    }
};

std::pair<float, std::vector<Nam>> find_nams(
    const QueryRandstrobeVector &query_randstrobes,
    const StrobemerIndex& index
);

std::vector<Nam> find_nams_rescue(
    const QueryRandstrobeVector &query_randstrobes,
    const StrobemerIndex& index,
    unsigned int rescue_cutoff
);

bool partial_search(
    std::vector<PartialSeed>& partial_queried,
    std::array<robin_hood::unordered_map<unsigned int, std::vector<MainHit>>, 2> &hits_per_ref,
    const StrobemerIndex& index,
    const QueryRandstrobe &qr,
    int &total_hits,
    int &nr_good_hits,
    uint aux_len,
    uint level
);

bool partial_search_rescue(
    std::vector<PartialSeed>& partial_queried,
    std::vector<RescueHit>& hits_fw,
    std::vector<RescueHit>& hits_rc,
    const StrobemerIndex& index,
    const QueryRandstrobe &qr,
    uint aux_len,
    uint level
);

std::ostream& operator<<(std::ostream& os, const Nam& nam);

#endif
