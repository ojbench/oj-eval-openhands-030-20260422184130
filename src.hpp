#ifndef SRC_HPP
#define SRC_HPP

#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include <cstddef>

typedef std::vector<std::vector<double>> IMAGE_T;

namespace nr_internal {
static inline void binarize(const IMAGE_T &img, std::vector<std::vector<unsigned char>> &bin, double th = 0.5) {
    const int n = (int)img.size();
    bin.assign(n, std::vector<unsigned char>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            bin[i][j] = img[i][j] > th ? 1 : 0;
        }
    }
}

static inline bool bbox(const std::vector<std::vector<unsigned char>> &b, int &r0, int &c0, int &r1, int &c1) {
    int n = (int)b.size();
    r0 = n; c0 = n; r1 = -1; c1 = -1;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (b[i][j]) {
                r0 = std::min(r0, i); c0 = std::min(c0, j);
                r1 = std::max(r1, i); c1 = std::max(c1, j);
            }
        }
    }
    if (r1 < r0 || c1 < c0) return false;
    return true;
}

struct HoleInfo {
    int count;
    double cx, cy; // centroid of merged holes (average)
};

static inline HoleInfo count_holes(const std::vector<std::vector<unsigned char>> &b) {
    int n = (int)b.size();
    // pad with 0 (background) around
    int H = n + 2, W = n + 2;
    std::vector<std::vector<unsigned char>> grid(H, std::vector<unsigned char>(W, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            grid[i+1][j+1] = b[i][j] ? 2 : 0; // 2 for foreground, 0 for background

    // mark outer background with 1 using BFS
    std::queue<std::pair<int,int>> q;
    q.push({0,0});
    grid[0][0] = 1;
    const int dr[4] = {1,-1,0,0};
    const int dc[4] = {0,0,1,-1};
    while(!q.empty()){
        auto [r,c] = q.front(); q.pop();
        for(int k=0;k<4;++k){
            int nr=r+dr[k], nc=c+dc[k];
            if(nr>=0 && nr<H && nc>=0 && nc<W && grid[nr][nc]==0){
                grid[nr][nc]=1; q.push({nr,nc});
            }
        }
    }
    // remaining 0s inside are holes
    int holes=0; long long sumr=0, sumc=0, cnt=0;
    for(int i=0;i<H;++i){
        for(int j=0;j<W;++j){
            if(grid[i][j]==0){
                // new hole component
                ++holes;
                std::queue<std::pair<int,int>> q2;
                q2.push({i,j}); grid[i][j]=3;
                while(!q2.empty()){
                    auto [r,c]=q2.front(); q2.pop();
                    sumr += (r-1); sumc += (c-1); ++cnt; // map back to original coords approx
                    for(int k=0;k<4;++k){
                        int nr=r+dr[k], nc=c+dc[k];
                        if(nr>=0 && nr<H && nc>=0 && nc<W && grid[nr][nc]==0){
                            grid[nr][nc]=3; q2.push({nr,nc});
                        }
                    }
                }
            }
        }
    }
    HoleInfo info; info.count = holes;
    if (cnt>0) { info.cy = (double)sumr / cnt; info.cx = (double)sumc / cnt; }
    else { info.cy = info.cx = -1.0; }
    return info;
}

static inline void projections_and_stats(const std::vector<std::vector<unsigned char>> &b,
                                         int r0,int c0,int r1,int c1,
                                         std::vector<int> &row, std::vector<int> &col,
                                         int &total, double &cx, double &cy) {
    int h = r1-r0+1, w = c1-c0+1;
    row.assign(h,0); col.assign(w,0); total=0;
    double sx=0.0, sy=0.0;
    for(int i=r0;i<=r1;++i){
        for(int j=c0;j<=c1;++j){
            if(b[i][j]){
                ++row[i-r0]; ++col[j-c0]; ++total;
                sy += (i - r0 + 0.5); sx += (j - c0 + 0.5);
            }
        }
    }
    if(total>0){ cy = sy/total; cx = sx/total; }
    else { cx = cy = 0.0; }
}

static inline double band_density_rows(const std::vector<int> &row, int a, int b){
    a = std::max(a,0); b = std::min(b,(int)row.size()); if(a>=b) return 0.0;
    long long s=0; for(int i=a;i<b;++i) s+=row[i];
    long long tot=0; for(int v:row) tot+=v;
    if(tot==0) return 0.0; return (double)s/(double)tot;
}
static inline double band_density_cols(const std::vector<int> &col, int a, int b){
    a = std::max(a,0); b = std::min(b,(int)col.size()); if(a>=b) return 0.0;
    long long s=0; for(int i=a;i<b;++i) s+=col[i];
    long long tot=0; for(int v:col) tot+=v;
    if(tot==0) return 0.0; return (double)s/(double)tot;
}

} // namespace

int judge(IMAGE_T &img) {
    using namespace nr_internal;
    if (img.empty() || img[0].empty()) return 0;
    std::vector<std::vector<unsigned char>> bin;
    binarize(img, bin, 0.5);
    int r0,c0,r1,c1; if(!bbox(bin,r0,c0,r1,c1)) return 0;

    int h = r1-r0+1, w = c1-c0+1;
    double wh_ratio = w>0 ? (double)w/(double)h : 0.0;

    HoleInfo hi = count_holes(bin);
    std::vector<int> row, col; int total=0; double cx=0, cy=0;
    projections_and_stats(bin, r0,c0,r1,c1, row, col, total, cx, cy);

    // quick path: distinguish by holes
    if (hi.count >= 2) return 8;
    if (hi.count == 1) {
        double hy = hi.cy - r0; // hole centroid relative to bbox
        double dx = std::abs((hi.cx - c0) - cx) / std::max(1.0, (double)w);
        double dy = std::abs(hy - cy) / std::max(1.0, (double)h);
        // 0: hole near center, more balanced shape
        if (dy < 0.15 && dx < 0.15) return 0;
        // 9 vs 6 by hole vertical position
        if (hy < cy) return 9; else return 6;
    }

    // No holes. Use shape heuristics
    // Densities
    double top = band_density_rows(row, 0, (int)std::round(0.25*row.size()));
    double bottom = band_density_rows(row, (int)std::round(0.75*row.size()), row.size());
    double mid = band_density_rows(row, (int)std::round(0.40*row.size()), (int)std::round(0.60*row.size()));
    double center_vert = band_density_cols(col, (int)std::round(0.4*col.size()), (int)std::round(0.6*col.size()));
    double leftd = band_density_cols(col, 0, (int)std::round(0.5*col.size()));
    double rightd = band_density_cols(col, (int)std::round(0.5*col.size()), col.size());

    // '1': narrow width and central vertical band
    if (wh_ratio < 0.55 || (center_vert > 0.65 && wh_ratio < 0.75 && top < 0.35)) {
        return 1;
    }

    // '7': strong top bar, light bottom, right-heavy
    if (top > 0.34 && bottom < 0.25 && rightd > leftd * 1.15) {
        return 7;
    }

    // '4': strong mid bar, bottom-left sparse, right heavy
    double bottom_left = 0.0, bottom_right = 0.0;
    {
        int rmid0 = (int)std::round(0.6*row.size());
        int c_mid = (int)std::round(0.5*col.size());
        bottom_left = band_density_rows(row, rmid0, row.size()) * band_density_cols(col, 0, c_mid);
        bottom_right = band_density_rows(row, rmid0, row.size()) * band_density_cols(col, c_mid, col.size());
    }
    if (mid > 0.28 && bottom_left < 0.06 && bottom_right > 0.10) {
        return 4;
    }

    // '3': left half very sparse and right heavy
    if (leftd < 0.20 && rightd > 0.55) {
        return 3;
    }

    // Distinguish 2 and 5 by horizontal centroid and bottom halves
    double cx_norm = cx / std::max(1.0, (double)w);
    double lower = band_density_rows(row, (int)std::round(0.5*row.size()), row.size());
    double lower_left = 0.0, lower_right = 0.0;
    {
        int rlow = (int)std::round(0.65*row.size());
        int c_mid = (int)std::round(0.5*col.size());
        // approximate by combining bands
        lower_left = band_density_rows(row, rlow, row.size()) * band_density_cols(col, 0, c_mid);
        lower_right = band_density_rows(row, rlow, row.size()) * band_density_cols(col, c_mid, col.size());
    }

    if (lower > 0.28) {
        if (lower_right > lower_left * 1.2 || cx_norm > 0.52) return 2; // more weight on right
        else return 5;
    }

    // Fallbacks by center of mass
    if (cx_norm > 0.55) return 2; // lean right
    if (cx_norm < 0.45) return 5; // lean left

    // Default guess
    return 2;
}

#endif // SRC_HPP
