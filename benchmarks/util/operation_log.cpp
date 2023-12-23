#include "operation_log.hpp"

#include "replay_tree.hpp"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <unordered_set>
#include <vector>

void operation_log::write(OperationLog const& logs, std::ostream& out) {
    out << logs.pushes.size() << ' ' << logs.pops.size() << '\n';
    auto pop_it = logs.pops.begin();
    auto push_it = logs.pushes.begin();
    while (pop_it != logs.pops.end() && push_it != logs.pushes.end()) {
        if (pop_it->tick < push_it->tick) {
            out << "d," << pop_it->tick << ',' << pop_it->ref_index << '\n';
            ++pop_it;
        } else {
            out << "i," << push_it->tick << ',' << push_it->key << ',' << push_it->index << '\n';
            ++push_it;
        }
    }
    while (pop_it != logs.pops.end()) {
        out << "d," << pop_it->tick << ',' << pop_it->ref_index << '\n';
        ++pop_it;
    }
    while (push_it != logs.pushes.end()) {
        out << "i," << push_it->tick << ',' << push_it->key << ',' << push_it->index << '\n';
        ++push_it;
    }
}
