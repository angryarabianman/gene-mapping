#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <cstring>
#include <cmath>
using namespace std;

const int BATCH_SIZE = 500000;
const int SEED_LENGTH = 16;
const int MIN_SEED_COUNT = 2;
const double SCORE_THRESHOLD_PCT = 0.85;

inline int get_id(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 0;
    }
}

string get_rc_seq(const string &seq) {
    string rc = "";
    rc.reserve(seq.size());
    for (int i = seq.size() - 1; i >= 0; i--) {
        char c = seq[i];
        if (c == 'A') rc += 'T';
        else if (c == 'T') rc += 'A';
        else if (c == 'C') rc += 'G';
        else if (c == 'G') rc += 'C';
        else rc += 'N';
    }
    return rc;
}

string get_rc_qual(const string &qual) {
    string rev = qual;
    reverse(rev.begin(), rev.end());
    return rev;
}

const int MAX_NODES = BATCH_SIZE * 5 * SEED_LENGTH + 1000;

struct SeedInfo {
    int read_index;
    int offset;
    bool is_rc;
};

struct AC {
    int next[MAX_NODES][4];
    int fail[MAX_NODES];
    vector<SeedInfo> outputs[MAX_NODES];
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
        for (int i = 0; i < 4; i++) {
            if (next[0][i]) q.push(next[0][i]);
        }
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

    double max_score = len;
    double current_score = 0.0;

    for (int i = 0; i < len; i++) {
        char g_char = genome[start_pos + i];
        char r_char = read[i];

        if (r_char == 'N' || g_char == 'N') { }
        else if (g_char == r_char) {
            current_score += 1.0;
        }
        else {
            int phred = qual[i] - 33;
            if (phred < 20) current_score -= 0.5;
            else current_score -= 3.0;
        }
    }

    return current_score >= (max_score * SCORE_THRESHOLD_PCT);
}

AC ac;
vector<int> coverage_diff;
long long total_reads = 0;
long long mapped_reads = 0;
long long unique_mapped = 0;
long long multi_mapped = 0;
double total_alignment_score = 0;
long long total_aligned_bases = 0;

struct ReadBatchData {
    string seq;
    string qual;
    string rc_seq;
    string rc_qual;
    bool mapped;
    int map_count;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string genome_file = "genome.fna";
    string reads_file = "reads.fastq";

    cout << "[1/4] Loading Genome..." << endl;
    ifstream f_gen(genome_file);
    if (!f_gen.is_open()) { cerr << "Error: genome.fna not found!" << endl; return 1; }

    string genome;
    string line;
    while(getline(f_gen, line)) {
        genome += line;
    }
    f_gen.close();
    int G_LEN = genome.size();
    coverage_diff.resize(G_LEN + 2, 0);

    cout << "      Genome Size: " << G_LEN << " bp" << endl;

    cout << "[2/4] Processing Reads (Batch Size: " << BATCH_SIZE << ")..." << endl;
    ifstream f_reads(reads_file);
    if (!f_reads.is_open()) { cerr << "Error: reads.fastq not found!" << endl; return 1; }

    // int cnt_lines = 0;
    // while(getline(f_reads, line)) {
    //     cnt_lines++;
    // }

    vector<ReadBatchData> batch;
    bool file_good = true;


    int cnt = 0;
    while (file_good) {
        // cnt++;
        // if (cnt > 10) break;
        batch.clear();
        ac.reset();

        for (int k = 0; k < BATCH_SIZE; k++) {
            string l1, seq, l3, qual;
            if (!getline(f_reads, l1)) { file_good = false; break; }
            getline(f_reads, seq);
            getline(f_reads, l3);
            getline(f_reads, qual);

            string rc_s = get_rc_seq(seq);
            string rc_q = get_rc_qual(qual);

            batch.push_back({seq, qual, rc_s, rc_q, false, 0});

            int idx = batch.size() - 1;

            if (seq.size() >= SEED_LENGTH)
                ac.insert_seed(seq.substr(0, SEED_LENGTH), idx, 0, false);
            if (seq.size() >= SEED_LENGTH * 2)
                ac.insert_seed(seq.substr(seq.size()/2, SEED_LENGTH), idx, seq.size()/2, false);

            if (rc_s.size() >= SEED_LENGTH)
                ac.insert_seed(rc_s.substr(0, SEED_LENGTH), idx, 0, true);
            if (rc_s.size() >= SEED_LENGTH * 2)
                ac.insert_seed(rc_s.substr(rc_s.size()/2, SEED_LENGTH), idx, rc_s.size()/2, true);
        }

        if (batch.empty()) break;
        ac.build_fail();

        int u = 0;
        for (int i = 0; i < G_LEN; i++) {
            char base = genome[i];
            int c = get_id(base);

            if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
                u = 0;
                continue;
            }

            u = ac.next[u][c];

            int temp = u;
            while (temp > 0) {
                for (const auto& match : ac.outputs[temp]) {
                    int r_idx = match.read_index;

                    if (batch[r_idx].map_count > 5) continue;
                    int read_start = (i - SEED_LENGTH + 1) - match.offset;
                    bool is_valid = false;
                    int r_len = batch[r_idx].seq.size();

                    if (match.is_rc) {
                        is_valid = verify_alignment(genome, read_start, batch[r_idx].rc_seq, batch[r_idx].rc_qual);
                    } else {
                        is_valid = verify_alignment(genome, read_start, batch[r_idx].seq, batch[r_idx].qual);
                    }

                    if (is_valid) {
                        if (batch[r_idx].map_count == 0) {
                            batch[r_idx].mapped = true;
                            if (read_start >= 0 && read_start + r_len < G_LEN) {
                                coverage_diff[read_start]++;
                                coverage_diff[read_start + r_len]--;
                            }
                        }
                        batch[r_idx].map_count++;
                    }
                }
                temp = ac.fail[temp];

                if (ac.outputs[temp].empty() && temp != 0) {
                    break;
                }
            }
        }

        for (const auto& r : batch) {
            total_reads++;
            if (r.mapped) {
                mapped_reads++;
                if (r.map_count == 1) unique_mapped++;
                else multi_mapped++;
            }
        }

        cout << "      Processed " << total_reads << " reads..." << endl;
    }

    cout << "[3/4] Calculating Coverage..." << endl;
    long long covered_bases = 0;
    long long current_depth = 0;
    for (int i = 0; i < G_LEN; i++) {
        current_depth += coverage_diff[i];
        if (current_depth > 0) covered_bases++;
    }

    cout << "\n================ REPORT ================" << endl;
    cout << "\n1. Statistics:" << endl;
    cout << "   - Total Reads:      " << total_reads << endl;
    cout << "   - Mapped Reads:     " << mapped_reads << " (" << (total_reads ? (double)mapped_reads/total_reads*100.0 : 0) << "%)" << endl;
    cout << "   - Unique Mapped:    " << unique_mapped << " (" << (mapped_reads ? (double)unique_mapped/mapped_reads*100.0 : 0) << "%)" << endl;
    cout << "   - Multi Mapped:     " << multi_mapped << " (" << (mapped_reads ? (double)multi_mapped/mapped_reads*100.0 : 0) << "%)" << endl;
    cout << "\n2. Quality:" << endl;
    cout << "   - Alignment Score Threshold: " << (SCORE_THRESHOLD_PCT*100) << "% (using Phred penalties)" << endl;
    cout << "\n3. Genome Coverage:" << endl;
    cout << "   - Covered Bases:    " << covered_bases << endl;
    cout << "   - Coverage %:       " << (G_LEN ? (double)covered_bases/G_LEN*100.0 : 0) << "%" << endl;
    cout << "========================================" << endl;

    return 0;
}