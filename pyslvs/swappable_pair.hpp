class SwappablePair {
    int first, second;

public:
    SwappablePair() = default;

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
