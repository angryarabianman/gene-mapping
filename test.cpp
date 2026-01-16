#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <cstring>
#include <cmath>

using namespace std;

const int SEED_LENGTH = 4;
const double SCORE_THRESHOLD_PCT = 0.85;

inline int get_id(char c) {
    switch (c) {
        case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return 0;
    }
}

string get_rc_seq(const string &seq) {
    string rc = "";
    for (int i = seq.size() - 1; i >= 0; i--) {
        char c = seq[i];
        if (c == 'A') rc += 'T'; else if (c == 'T') rc += 'A';
        else if (c == 'C') rc += 'G'; else if (c == 'G') rc += 'C'; else rc += 'N';
    }
    return rc;
}

string get_rc_qual(const string &qual) {
    string rev = qual;
    reverse(rev.begin(), rev.end());
    return rev;
}

struct SeedInfo { int read_index; int offset; bool is_rc; };

struct AC {
    static const int MAX_TEST_NODES = 10000;
    int next[MAX_TEST_NODES][4];
    int fail[MAX_TEST_NODES];
    vector<SeedInfo> outputs[MAX_TEST_NODES];
    int nodes_cnt;

    AC() { reset(); }

    void reset() {
        nodes_cnt = 1;
        memset(next[0], 0, sizeof(next[0]));
        fail[0] = 0;
        outputs[0].clear();
    }

    int new_node() {
        int u = nodes_cnt++;
        memset(next[u], 0, sizeof(next[u]));
        fail[u] = 0;
        outputs[u].clear();
        return u;
    }

    void insert_seed(const string &s, int read_idx, int offset, bool rc) {
        int u = 0;
        for (char c : s) {
            if (c == 'N') return;
            int v = get_id(c);
            if (!next[u][v]) next[u][v] = new_node();
            u = next[u][v];
        }
        outputs[u].push_back({read_idx, offset, rc});
    }

    void build_fail() {
        queue<int> q;
        for (int i = 0; i < 4; i++) if (next[0][i]) q.push(next[0][i]);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int i = 0; i < 4; i++) {
                if (next[u][i]) {
                    fail[next[u][i]] = next[fail[u]][i];
                    q.push(next[u][i]);
                } else {
                    next[u][i] = next[fail[u]][i];
                }
            }
        }
    }
};

bool verify_alignment(const string &genome, int start_pos, 
                      const string &read, const string &qual) {
    int len = read.size();
    if (start_pos < 0 || start_pos + len > (int)genome.size()) return false;

    double max_score = (double)len; 
    double current_score = 0.0;

    for (int i = 0; i < len; i++) {
        char g_char = genome[start_pos + i];
        char r_char = read[i];
        
        if (r_char == 'N' || g_char == 'N') { }
        else if (g_char == r_char) { current_score += 1.0; } 
        else {
            int phred = qual[i] - 33;
            if (phred < 20) current_score -= 0.5;
            else current_score -= 3.0;
        }
    }
    return current_score >= (max_score * SCORE_THRESHOLD_PCT);
}

void run_test_case(string test_name, string genome, string read, string qual, bool expect_map) {
    cout << "TEST: " << test_name << "... ";
    
    AC ac;
    if (read.size() >= SEED_LENGTH) 
        ac.insert_seed(read.substr(0, SEED_LENGTH), 0, 0, false);
    if (read.size() >= SEED_LENGTH * 2) 
        ac.insert_seed(read.substr(read.size()/2, SEED_LENGTH), 0, read.size()/2, false);

    string rc_s = get_rc_seq(read);
    string rc_q = get_rc_qual(qual);
    if (rc_s.size() >= SEED_LENGTH) 
        ac.insert_seed(rc_s.substr(0, SEED_LENGTH), 0, 0, true);
    if (rc_s.size() >= SEED_LENGTH * 2) 
        ac.insert_seed(rc_s.substr(rc_s.size()/2, SEED_LENGTH), 0, rc_s.size()/2, true);

    ac.build_fail();

    bool mapped = false;
    int u = 0;
    for (int i = 0; i < genome.size(); i++) {
        char base = genome[i];
        if (base == 'N') { u = 0; continue; }
        u = ac.next[u][get_id(base)];
        
        int temp = u;
        while (temp > 0) {
            for (const auto& match : ac.outputs[temp]) {
                int read_start = (i - SEED_LENGTH + 1) - match.offset;
                
                bool result = false;
                if (match.is_rc) result = verify_alignment(genome, read_start, rc_s, rc_q);
                else result = verify_alignment(genome, read_start, read, qual);

                if (result) mapped = true;
            }
            temp = ac.fail[temp];
            if (mapped) break;
        }
        if (mapped) break;
    }

    if (mapped == expect_map) {
        cout << "PASSED" << endl;
    } else {
        cout << "FAILED!" << endl;
        cout << "   Genome: " << genome << endl;
        cout << "   Read:   " << read << endl;
        cout << "   Expect: " << (expect_map ? "Map" : "No Map") << endl;
        exit(1);
    }
}

int main() {
    cout << "=== Running Unit Tests ===" << endl;

    string high_qual = "IIIIIIIIIIIIIIII";
    string low_qual  = "################";

    run_test_case("Exact Match Forward",
                  "CCCAAAACCCC",
                  "AAAA",
                  "IIII",
                  true);

    run_test_case("Reverse Complement Match",
                  "CCCTTTTCCCC",
                  "AAAA",
                  "IIII",
                  true);

    run_test_case("Low Quality Mismatch",
                  "CCCAAAAAAAAAAAACCCC",
                  "AAAAAAAAAAAT",
                  "IIIIIIIIIII#",
                  true);

    run_test_case("High Quality Mismatch",
                  "CCCAAAACCCC",
                  "AAAT",
                  "IIII",
                  false);

    run_test_case("Seed Rescue",
                  "GCAAAAAAAAAAAACCCCGT",
                  "AAAAAAAAAAATCCCC",
                  "IIIIIIIIIII#IIII",
                  true);

    cout << "=== All Tests Passed ===" << endl;
    return 0;
}