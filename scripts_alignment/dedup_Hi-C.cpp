#include <cmath>
#include <ctime>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <sys/resource.h>
using namespace std;

int num_shift_R1, num_shift_R2, num_cell, num_chrom;

map<string,int> chrom_a2i, cell_a2i;
vector<string> chrom_i2a, cell_i2a;
vector<vector<pair<int,int>>> cc2rp;

int encode_chrom(char *chrom, char strand) {
	string c = chrom;
	auto it = chrom_a2i.find(c);
	if (it == chrom_a2i.end()) {
		it = chrom_a2i.insert(make_pair(c, chrom_a2i.size())).first;
		chrom_i2a.push_back(c);
		assert(chrom_a2i.size() <= num_chrom);
	}
	return (it->second << 1) + (strand == '+'? 0: 1);
}

void decode_chrom(int x, string &c, char &s) {
	s = x&1? '-': '+';
	c = chrom_i2a[x>>1];
}

int encode_cell(char *cell) {
	string c = cell;
	auto it = cell_a2i.find(c);
	if (it == cell_a2i.end()) {
		it = cell_a2i.insert(make_pair(c, cell_a2i.size())).first;
		cell_i2a.push_back(c);
		assert(cell_a2i.size() <= num_cell);
	}
	return it->second;
}

string decode_cell(int x) {
	return cell_i2a[x];
}

int encode_cell_chrom(int cell, int c1, int c2) {
	return (cell * num_chrom * 2 + c1) * num_chrom * 2 + c2;
}

bool decode_cell_chrom(int key, int &cell, int &c1, int &c2) {
	c2 = key % (num_chrom * 2); key /= (num_chrom * 2);
	c1 = key % (num_chrom * 2); key /= (num_chrom * 2);
	cell = key;
	return cell < cell_i2a.size() && c1 < chrom_i2a.size()*2 && c2 < chrom_i2a.size()*2;
}

void insert_read_pair(int cell, int c1, int c2, int pos1, int pos2) {
	int key = encode_cell_chrom(cell, c1, c2);
	cc2rp[key].push_back(make_pair(pos1, pos2));
}

class UnionFindSet {
private:
public:
	vector<int> f;
	UnionFindSet(int n) {f.resize(n, -1);}
	int find(int x) { return f[x] < 0? x: f[x] = find(f[x]); }
	bool merge(int x, int y) {
		if ((x = find(x)) == (y = find(y))) return false;
		if (f[x] > f[y]) swap(x, y);
		f[x] += f[y]; f[y] = x;
		return true;
	}
};

bool keep_first(const pair<int,int> &p, const pair<int,int> &q) {
	int dp = abs(p.first-p.second), dq = abs(q.first-q.second);
	return dp != dq? dp > dq: min(p.first, p.second) <= min(q.first, q.second);
}

double get_duration(chrono::time_point<chrono::high_resolution_clock> _t_begin) {
// return time in second
	return chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now()-_t_begin).count() * 1e-3;
}

int get_mem_usage() {
// in KB
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	return usage.ru_maxrss;
}

vector<int> dedup(vector<pair<int,int>> &pairs) {
	int n = pairs.size();
	UnionFindSet ufs(n);
	sort(pairs.begin(), pairs.end());
	set<pair<int,int>> pos_id;
	for (int l = 0, i = 0; i < n; i++) {
		for (; l < i && pairs[l].first + num_shift_R1 < pairs[i].first; l++) pos_id.erase(make_pair(pairs[l].second, l));
		auto it_l = pos_id.lower_bound(make_pair(pairs[i].second-num_shift_R2, 0));
		auto it_r = pos_id.upper_bound(make_pair(pairs[i].second+num_shift_R2, n));
		for (auto it = it_l; it != it_r; it++) ufs.merge(i, it->second);
		pos_id.insert(make_pair(pairs[i].second, i));
	}
	for (int i = 0, j; i < n; i++) if ((j = ufs.find(i)) != i && keep_first(pairs[i], pairs[j])) {
		ufs.f[i] = ufs.f[j]; ufs.f[j] = i;
	}
	return ufs.f;
}

int main(int argc, char** argv) {
	num_shift_R1 = atoi(argv[1]);
	num_shift_R2 = atoi(argv[2]);
	num_cell = atoi(argv[3]);
	num_chrom = atoi(argv[4]);

	cc2rp.resize(num_cell * num_chrom * num_chrom * 4);
	int mem_usage_init = get_mem_usage();
	fprintf(stderr, "mem usage %.3fGB\n", mem_usage_init * 1e-6);

	{
	auto _t_begin = chrono::high_resolution_clock::now();
	char chrom1[32], chrom2[32], cell_raw[16], strand1, strand2, flipped;
	int pos1, pos2;
	int cnt_rp = 0;
	while(scanf(
		"%s\t%d\t%s\t%d\t%c\t%c\t%s\t%c\n",
		chrom1, &pos1, chrom2, &pos2, &strand1, &strand2, cell_raw, &flipped) != EOF
	) {
		int c1 = encode_chrom(chrom1, strand1), c2 = encode_chrom(chrom2, strand2);
		int cell = encode_cell(cell_raw);
		insert_read_pair(cell, c1, c2, pos1, pos2);
		cnt_rp++;
		if (cnt_rp % 10000 == 0) {
			double time_passed = get_duration(_t_begin);
			int mem_usage = get_mem_usage();
			fprintf(stderr, (string("\r") + string(110, ' ')).c_str());
			fprintf(
				stderr, "\rtime passed %.1fmin\tloaded %d lines\t%.1f lines/s \tmem %.3fGB\t%dB/line",
				time_passed/60, cnt_rp, cnt_rp / time_passed, mem_usage * 1e-6, (int)((mem_usage - mem_usage_init) * 1e3 / cnt_rp));
		}
	}
	fprintf(stderr, "\n");
	}

	{
	auto _t_begin = chrono::high_resolution_clock::now();
	int num_items = cc2rp.size(), num_processed = 0;
	for (int key = 0; key < num_items; key++) {
		num_processed++;
		if (num_processed % 10000 == 0 || key+1 == num_items) {
			double time_passed = get_duration(_t_begin);
			double num_item_per_sec = num_processed / time_passed;
			double time_left = (num_items - num_processed) / num_item_per_sec;
			double time_total = num_items / num_item_per_sec;
			fprintf(stderr, (string("\r") + string(160, ' ')).c_str());
			fprintf(
				stderr, "\rtime passed %.1fmin\tprocessing %d/%d\t%.1f items/s\ttime left %.1fmin\ttime total %.1fmin",
				time_passed/60, num_processed, num_items, num_item_per_sec, time_left/60, time_total/60);
		}

		int celli, c1i, c2i;
		string cell, c1, c2; char s1, s2;
		if (!decode_cell_chrom(key, celli, c1i, c2i)) {assert(cc2rp[key].empty()); continue;}
		auto &rp = cc2rp[key];
		if (rp.empty()) continue;
		cell = decode_cell(celli);
		decode_chrom(c1i, c1, s1);
		decode_chrom(c2i, c2, s2);
		vector<int> ufs_data = dedup(rp);
		for (int i = 0; i < rp.size(); i++) {
			const auto &p = rp[i]; int f = ufs_data[i];
			printf(
				".\t%s\t%d\t%s\t%d\t%c\t%c\t%s\t%s\t%d\n",
				c1.c_str(), p.first, c2.c_str(), p.second, s1, s2,
				f<0? "..": "DD", cell.c_str(), f<0? -f: 0);
		}
	}
	fprintf(stderr, "\n");
	}
	return 0;
}
