#include "overlap.hpp"
#include <vector>
#include <unordered_map>
#include <tuple>
#include <algorithm>
#include <cassert>
using namespace std;

struct State {
    int len, link;
    int first_pos; // posizione di fine della prima occorrenza
    unordered_map<unsigned int, int> next;
};

struct SuffixAutomaton {
    vector<State> st;
    int last;
};

static SuffixAutomaton build_sa(const vector<unsigned int>& A) {
    SuffixAutomaton sa;
    sa.st.reserve(2 * A.size());
    State initial;
    initial.len = 0;
    initial.link = -1;
    initial.first_pos = 0;
    sa.st.push_back(initial);
    sa.last = 0;
    for (int i = 0; i < A.size(); i++) {
        unsigned int c = A[i];
        int cur = sa.st.size();
        State curState;
        curState.len = sa.st[sa.last].len + 1;
        curState.first_pos = i;
        curState.link = 0;
        sa.st.push_back(curState);
        int p = sa.last;
        for (; p != -1 && sa.st[p].next.find(c) == sa.st[p].next.end(); p = sa.st[p].link) {
            sa.st[p].next[c] = cur;
        }
        if (p == -1) {
            sa.st[cur].link = 0;
        } else {
            int q = sa.st[p].next[c];
            if (sa.st[p].len + 1 == sa.st[q].len) {
                sa.st[cur].link = q;
            } else {
                int clone = sa.st.size();
                State cloneState;
                cloneState.len = sa.st[p].len + 1;
                cloneState.next = sa.st[q].next; // copia della mappa
                cloneState.link = sa.st[q].link;
                cloneState.first_pos = sa.st[q].first_pos;
                sa.st.push_back(cloneState);
                for (; p != -1 && sa.st[p].next[c] == q; p = sa.st[p].link) {
                    sa.st[p].next[c] = clone;
                }
                sa.st[q].link = sa.st[cur].link = clone;
            }
        }
        sa.last = cur;
    }
    return sa;
}

tuple<int, int, int> longest_common_substring_suffix_automaton(const vector<unsigned int>& A, const vector<unsigned int>& B) {
    SuffixAutomaton sa = build_sa(A);
    int v = 0, l = 0;
    int best = 0, bestposB = 0;
    int best_state = 0;
    for (int i = 0; i < B.size(); i++) {
        unsigned int c = B[i];
        if (sa.st[v].next.find(c) != sa.st[v].next.end()) {
            v = sa.st[v].next[c];
            l++;
        } else {
            while (v != -1 && sa.st[v].next.find(c) == sa.st[v].next.end())
                v = sa.st[v].link;
            if (v == -1) {
                v = 0;
                l = 0;
            } else {
                l = sa.st[v].len + 1;
                v = sa.st[v].next[c];
            }
        }
        if (l > best) {
            best = l;
            bestposB = i;
            best_state = v;
        }
    }
    int startB = bestposB - best + 1;
    int startA = sa.st[best_state].first_pos - best + 1;
    return make_tuple(best, startA, startB);
}
