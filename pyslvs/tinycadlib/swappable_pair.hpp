#ifndef SWAPPABLE_PAIR_HPP
#define SWAPPABLE_PAIR_HPP

struct SwappablePair {
    int first, second;

    inline auto operator==(const SwappablePair &rhs) const -> bool {
        return (first == rhs.first && second == rhs.second) ||
               (first == rhs.second && second == rhs.first);
    }

    inline auto operator!=(const SwappablePair &rhs) const -> bool {
        return !this->operator==(rhs);
    }

    inline auto operator<(const SwappablePair &rhs) const -> bool {
        return this->operator!=(rhs) && (first < rhs.first || second < rhs.second);
    }
};

#endif //SWAPPABLE_PAIR_HPP
