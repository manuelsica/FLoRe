#include "overlap.hpp"
#include <algorithm>

SuffixAutomaton build_suffix_automaton(const std::vector<unsigned int>& A) {
    SuffixAutomaton sa;
    sa.st.reserve(2 * A.size());
    sa.last = 0;

    // Stato iniziale
    State initial;
    initial.len = 0;
    initial.link = -1;
    initial.first_pos = 0;
    sa.st.push_back(initial);

    // Costruzione
    for (int i = 0; i < (int)A.size(); i++){
        unsigned int c = A[i];
        int cur = (int)sa.st.size();
        State curState;
        curState.len = sa.st[sa.last].len + 1;
        curState.link = 0;
        curState.first_pos = i;
        sa.st.push_back(curState);

        int p = sa.last;
        while (p != -1 && sa.st[p].next.find(c) == sa.st[p].next.end()){
            sa.st[p].next[c] = cur;
            p = sa.st[p].link;
        }
        if (p == -1){
            sa.st[cur].link = 0;
        } else {
            int q = sa.st[p].next[c];
            if (sa.st[p].len + 1 == sa.st[q].len){
                sa.st[cur].link = q;
            } else {
                int clone = (int)sa.st.size();
                State cloneState;
                cloneState.len = sa.st[p].len + 1;
                cloneState.next = sa.st[q].next;
                cloneState.link = sa.st[q].link;
                cloneState.first_pos = sa.st[q].first_pos;
                sa.st.push_back(cloneState);
                while (p != -1 && sa.st[p].next[c] == q){
                    sa.st[p].next[c] = clone;
                    p = sa.st[p].link;
                }
                sa.st[q].link = clone;
                sa.st[cur].link = clone;
            }
        }
        sa.last = cur;
    }
    return sa;
}

std::tuple<int,int,int> match_suffix_automaton(const SuffixAutomaton &sa,
                                               const std::vector<unsigned int> &B)
{
    int v = 0;
    int l = 0;
    int best = 0;
    int bestposB = 0;
    int best_state = 0;

    for(int i = 0; i < (int)B.size(); i++){
        unsigned int c = B[i];
        if(sa.st[v].next.find(c) != sa.st[v].next.end()){
            v = sa.st[v].next.at(c);
            l++;
        } else {
            while(v != -1 && sa.st[v].next.find(c) == sa.st[v].next.end()){
                v = sa.st[v].link;
            }
            if(v == -1){
                v = 0;
                l = 0;
            } else {
                l = sa.st[v].len + 1;
                v = sa.st[v].next.at(c);
            }
        }
        if(l > best){
            best = l;
            bestposB = i;
            best_state = v;
        }
    }
    int startB = bestposB - best + 1;
    int startA = sa.st[best_state].first_pos - best + 1;
    return std::make_tuple(best, startA, startB);
}

std::tuple<int,int,int> longest_common_substring_suffix_automaton(
    const std::vector<unsigned int>& A,
    const std::vector<unsigned int>& B)
{
    SuffixAutomaton sa = build_suffix_automaton(A);
    return match_suffix_automaton(sa, B);
}
