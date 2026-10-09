#pragma once

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cctype>
#include <charconv>
#include <cstddef>
#include <filesystem>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Reads DIMACS shortest path files ("p sp n m", arcs "a u v w") and METIS files ("n m [fmt [ncon]]" followed by one
// adjacency line per node, edge weights 1 if the format has none).
struct Graph {
    using weight_type = long long;
    struct Edge {
        std::size_t target;
        weight_type weight;
    };
    std::vector<std::size_t> nodes;
    std::vector<Edge> edges;

    Graph() = default;

    explicit Graph(std::filesystem::path const& graph_file) {
        int fd = open(graph_file.c_str(), O_RDONLY);
        if (fd == -1) {
            throw std::runtime_error{"Could not open file"};
        }
        struct stat sb{};
        if (fstat(fd, &sb) == -1) {
            close(fd);
            throw std::runtime_error{"Could not get file size"};
        }
        auto size = static_cast<std::size_t>(sb.st_size);
        if (size == 0) {
            close(fd);
            throw std::runtime_error{"Empty file"};
        }
        auto* addr = mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0);
        close(fd);
        if (addr == MAP_FAILED) {
            throw std::runtime_error{"mmap failed"};
        }
        madvise(addr, size, MADV_SEQUENTIAL);
        const auto* begin = static_cast<char const*>(addr);
        try {
            const auto* first = begin;
            while (first != begin + size && std::isspace(static_cast<unsigned char>(*first)) != 0) {
                ++first;
            }
            if (first != begin + size && (*first == 'c' || *first == 'p')) {
                read_dimacs(first, begin + size);
            } else {
                read_metis(first, begin + size);
            }
        } catch (...) {
            munmap(addr, size);
            throw;
        }
        munmap(addr, size);
    }

    [[nodiscard]] std::size_t num_nodes() const noexcept {
        return nodes.empty() ? 0 : nodes.size() - 1;
    }

    [[nodiscard]] std::size_t num_edges() const noexcept {
        return edges.size();
    }

    [[nodiscard]] std::size_t degree(std::size_t node) const noexcept {
        return nodes[node + 1] - nodes[node];
    }

    [[nodiscard]] std::size_t source(std::string const& spec) const {
        if (spec == "max-degree") {
            std::size_t best = 0;
            for (std::size_t i = 1; i < num_nodes(); ++i) {
                if (degree(i) > degree(best)) {
                    best = i;
                }
            }
            return best;
        }
        std::size_t node = 0;
        auto res = std::from_chars(spec.data(), spec.data() + spec.size(), node);
        if (res.ec != std::errc{} || res.ptr != spec.data() + spec.size()) {
            throw std::invalid_argument{"Source must be a node index or max-degree"};
        }
        if (node >= num_nodes()) {
            throw std::out_of_range{"Source node out of range"};
        }
        return node;
    }

   private:
    template <typename T>
    static char const* parse_number(char const* it, char const* end, T& value, char const* what) {
        auto res = std::from_chars(it, end, value);
        if (res.ec != std::errc{}) {
            throw std::runtime_error{std::string{"Failed to parse "} + what};
        }
        return res.ptr;
    }

    static char const* skip_space(char const* it, char const* end) {
        while (it != end && std::isspace(static_cast<unsigned char>(*it)) != 0) {
            ++it;
        }
        return it;
    }

    static char const* skip_blank(char const* it, char const* end) {
        while (it != end && (*it == ' ' || *it == '\t' || *it == '\r')) {
            ++it;
        }
        return it;
    }

    static char const* skip_line(char const* it, char const* end) {
        while (it != end && *it++ != '\n') {
        }
        return it;
    }

    void read_dimacs(char const* it, char const* end) {
        std::size_t num_nodes = 0;
        std::size_t num_edges = 0;
        bool header = false;
        std::vector<std::pair<std::size_t, Edge>> edge_list;
        while ((it = skip_space(it, end)) != end) {
            if (*it == 'c') {
                it = skip_line(it, end);
            } else if (*it == 'p') {
                if (header) {
                    throw std::runtime_error{"Duplicate problem line"};
                }
                it = skip_blank(it + 1, end);
                while (it != end && std::isspace(static_cast<unsigned char>(*it)) == 0) {
                    ++it;
                }
                it = parse_number(skip_blank(it, end), end, num_nodes, "number of nodes");
                it = parse_number(skip_blank(it, end), end, num_edges, "number of edges");
                header = true;
                nodes.assign(num_nodes + 1, 0);
                edge_list.reserve(num_edges);
            } else if (*it == 'a' && header) {
                std::pair<std::size_t, Edge> edge;
                it = parse_number(skip_blank(it + 1, end), end, edge.first, "edge source");
                it = parse_number(skip_blank(it, end), end, edge.second.target, "edge target");
                it = parse_number(skip_blank(it, end), end, edge.second.weight, "edge weight");
                if (edge.first == 0 || edge.first > num_nodes || edge.second.target == 0 ||
                    edge.second.target > num_nodes) {
                    throw std::runtime_error{"Edge endpoint out of range"};
                }
                --edge.first;
                --edge.second.target;
                ++nodes[edge.first + 1];
                edge_list.push_back(edge);
            } else {
                throw std::runtime_error{"Invalid line"};
            }
        }
        if (!header) {
            throw std::runtime_error{"Missing problem line"};
        }
        if (edge_list.size() != num_edges) {
            throw std::runtime_error{"Number of edges does not match the problem line"};
        }
        std::exclusive_scan(nodes.begin() + 1, nodes.end(), nodes.begin() + 1, std::size_t{0});
        edges.resize(edge_list.size());
        for (auto const& edge : edge_list) {
            edges[nodes[edge.first + 1]++] = edge.second;
        }
    }

    void read_metis(char const* it, char const* end) {
        auto skip_comments = [&] {
            while (it != end && *it == '%') {
                it = skip_line(it, end);
            }
        };
        skip_comments();
        std::size_t num_nodes = 0;
        std::size_t num_edges = 0;
        it = parse_number(skip_blank(it, end), end, num_nodes, "number of nodes");
        it = parse_number(skip_blank(it, end), end, num_edges, "number of edges");
        std::string fmt;
        it = skip_blank(it, end);
        while (it != end && std::isdigit(static_cast<unsigned char>(*it)) != 0) {
            fmt += *it++;
        }
        if (fmt.size() > 3) {
            throw std::runtime_error{"Invalid format field"};
        }
        fmt.insert(0, 3 - fmt.size(), '0');
        bool vertex_sizes = fmt[0] == '1';
        bool vertex_weights = fmt[1] == '1';
        bool edge_weights = fmt[2] == '1';
        std::size_t ncon = 1;
        it = skip_blank(it, end);
        if (it != end && *it != '\n') {
            it = parse_number(it, end, ncon, "number of constraints");
        }
        it = skip_line(it, end);

        nodes.assign(num_nodes + 1, 0);
        edges.reserve(2 * num_edges);
        for (std::size_t u = 0; u < num_nodes; ++u) {
            skip_comments();
            std::size_t skip = (vertex_sizes ? 1 : 0) + (vertex_weights ? ncon : 0);
            while ((it = skip_blank(it, end)) != end && *it != '\n') {
                long long value = 0;
                it = parse_number(it, end, value, "adjacency entry");
                if (skip > 0) {
                    --skip;
                    continue;
                }
                if (value <= 0 || static_cast<std::size_t>(value) > num_nodes) {
                    throw std::runtime_error{"Neighbor out of range"};
                }
                Edge edge{static_cast<std::size_t>(value - 1), 1};
                if (edge_weights) {
                    it = parse_number(skip_blank(it, end), end, edge.weight, "edge weight");
                }
                edges.push_back(edge);
            }
            if (it != end) {
                ++it;
            }
            nodes[u + 1] = edges.size();
        }
        if (edges.size() != 2 * num_edges) {
            throw std::runtime_error{"Number of adjacency entries does not match twice the number of edges"};
        }
    }
};
