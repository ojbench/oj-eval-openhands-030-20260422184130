#ifndef SRC_HPP
#define SRC_HPP

#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <utility>

namespace nr_internal {
static inline void binarize(const std::vector<std::vector<double> > &img, std::vector<std::vector<unsigned char> > &bin, double th = 0.5) {
    const int n = (int)img.size();
    bin.assign(n, std::vector<unsigned char>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            bin[i][j] = img[i][j] > th ? 1 : 0;
        }
    }
}

static inline bool bbox(const std::vector<std::vector<unsigned char> > &b, int &r0, int &c0, int &r1, int &c1) {
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

static inline std::vector<std::vector<unsigned char> > dilate1(const std::vector<std::vector<unsigned char> > &b){
    int n = (int)b.size();
    std::vector<std::vector<unsigned char> > out(n, std::vector<unsigned char>(n,0));
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j){
            if(!b[i][j]){
                bool on=false;
                for(int di=-1; di<=1 && !on; ++di){
                    for(int dj=-1; dj<=1 && !on; ++dj){
                        int ni=i+di, nj=j+dj;
                        if(ni>=0 && ni<n && nj>=0 && nj<n && b[ni][nj]) on=true;
                    }
                }
                out[i][j] = on ? 1 : 0;
            } else out[i][j]=1;
        }
    }
    return out;
}

static inline HoleInfo count_holes(const std::vector<std::vector<unsigned char> > &b) {
    int n = (int)b.size();
    int H = n + 2, W = n + 2;
    std::vector<std::vector<unsigned char> > grid(H, std::vector<unsigned char>(W, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            grid[i+1][j+1] = b[i][j] ? 2 : 0;

    std::queue<std::pair<int,int> > q;
    q.push(std::make_pair(0,0));
    grid[0][0] = 1;
    const int dr[4] = {1,-1,0,0};
    const int dc[4] = {0,0,1,-1};
    while(!q.empty()){
        std::pair<int,int> p = q.front(); q.pop();
        int r = p.first, c = p.second;
        for(int k=0;k<4;++k){
            int nr=r+dr[k], nc=c+dc[k];
            if(nr>=0 && nr<H && nc>=0 && nc<W && grid[nr][nc]==0){
                grid[nr][nc]=1; q.push(std::make_pair(nr,nc));
            }
        }
    }
    int holes=0; double sumr=0.0, sumc=0.0; int cnt=0;
    const int min_area = 5; // filter tiny cavities
    for(int i=0;i<H;++i){
        for(int j=0;j<W;++j){
            if(grid[i][j]==0){
                // new hole component
                int comp_cnt=0; double comp_sr=0.0, comp_sc=0.0;
                std::queue<std::pair<int,int> > q2;
                q2.push(std::make_pair(i,j)); grid[i][j]=3;
                while(!q2.empty()){
                    std::pair<int,int> p2 = q2.front(); q2.pop();
                    int r = p2.first, c = p2.second;
                    comp_sr += (r-1); comp_sc += (c-1); ++comp_cnt;
                    for(int k=0;k<4;++k){
                        int nr=r+dr[k], nc=c+dc[k];
                        if(nr>=0 && nr<H && nc>=0 && nc<W && grid[nr][nc]==0){
                            grid[nr][nc]=3; q2.push(std::make_pair(nr,nc));
                        }
                    }
                }
                if (comp_cnt >= min_area) {
                    ++holes;
                    sumr += comp_sr; sumc += comp_sc; cnt += comp_cnt;
                }
            }
        }
    }
    HoleInfo info; info.count = holes;
    if (cnt>0) { info.cy = sumr / cnt; info.cx = sumc / cnt; }
    else { info.cy = info.cx = -1.0; }
    return info;
}

static inline void projections_and_stats(const std::vector<std::vector<unsigned char> > &b,
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
    if (a < 0) a = 0; if (b > (int)row.size()) b = (int)row.size(); if(a>=b) return 0.0;
    double s=0.0; for(int i=a;i<b;++i) s+=row[i];
    double tot=0.0; for(size_t i=0;i<row.size();++i) tot+=row[i];
    if(tot==0.0) return 0.0; return s/tot;
}
static inline double band_density_cols(const std::vector<int> &col, int a, int b){
    if (a < 0) a = 0; if (b > (int)col.size()) b = (int)col.size(); if(a>=b) return 0.0;
    double s=0.0; for(int i=a;i<b;++i) s+=col[i];
    double tot=0.0; for(size_t i=0;i<col.size();++i) tot+=col[i];
    if(tot==0.0) return 0.0; return s/tot;
}

} // namespace

int judge(std::vector<std::vector<double> > &img) {
    using namespace nr_internal;
    if (img.empty() || img[0].empty()) return 0;
    std::vector<std::vector<unsigned char> > bin;
    binarize(img, bin, 0.5);
    int r0,c0,r1,c1; if(!bbox(bin,r0,c0,r1,c1)) return 0;

    int h = r1-r0+1, w = c1-c0+1;
    double wh_ratio = w>0 ? (double)w/(double)h : 0.0;

    HoleInfo hi = count_holes(bin);
    std::vector<int> row, col; int total=0; double cx=0, cy=0;
    projections_and_stats(bin, r0,c0,r1,c1, row, col, total, cx, cy);

    if (hi.count >= 2) return 8;
    if (hi.count == 1) {
        double hy = hi.cy - r0;
        double dx = std::fabs((hi.cx - c0) - cx) / (w>1 ? (double)w : 1.0);
        double dy = std::fabs(hy - cy) / (h>1 ? (double)h : 1.0);
        if (dy < 0.10 && dx < 0.10 && wh_ratio > 0.7 && wh_ratio < 1.3) return 0;
        if (hy < cy) return 9; else return 6;
    }

    int rs = (int)row.size(); int cs = (int)col.size();
    int t_end = (int)(0.25*rs + 0.5);
    int b_start = (int)(0.75*rs + 0.5);
    int m_start = (int)(0.40*rs + 0.5);
    int m_end = (int)(0.60*rs + 0.5);
    int cv_start = (int)(0.40*cs + 0.5);
    int cv_end = (int)(0.60*cs + 0.5);
    int c_mid = (int)(0.50*cs + 0.5);

    double top = band_density_rows(row, 0, t_end);
    double bottom = band_density_rows(row, b_start, rs);
    double mid = band_density_rows(row, m_start, m_end);
    double center_vert = band_density_cols(col, cv_start, cv_end);
    double leftd = band_density_cols(col, 0, c_mid);
    double rightd = band_density_cols(col, c_mid, cs);

    if (wh_ratio < 0.55 || (center_vert > 0.65 && wh_ratio < 0.75 && top < 0.35)) {
        return 1;
    }

    if (top > 0.34 && bottom < 0.25 && rightd > leftd * 1.15) {
        return 7;
    }

    double bottom_left = 0.0, bottom_right = 0.0;
    {
        int rmid0 = (int)(0.6*rs + 0.5);
        bottom_left = band_density_rows(row, rmid0, rs) * band_density_cols(col, 0, c_mid);
        bottom_right = band_density_rows(row, rmid0, rs) * band_density_cols(col, c_mid, cs);
    }
    if (mid > 0.35 && bottom < 0.18 && top > 0.12) {
        return 4;
    }

    if ((leftd < 0.30 && rightd > 0.45) && (top > bottom + 0.05)) {
        return 3;
    }

    double cx_norm = cx / (w>1 ? (double)w : 1.0);
    double lower = band_density_rows(row, (int)(0.5*rs + 0.5), rs);
    double lower_left = 0.0, lower_right = 0.0;
    {
        int rlow = (int)(0.65*rs + 0.5);
        lower_left = band_density_rows(row, rlow, rs) * band_density_cols(col, 0, c_mid);
        lower_right = band_density_rows(row, rlow, rs) * band_density_cols(col, c_mid, cs);
    }

    if (lower > 0.28) {
        if (lower_right > lower_left * 1.2 || cx_norm > 0.52) return 2;
        else return 5;
    }

    if (cx_norm > 0.55) return 2;
    if (cx_norm < 0.45) return 5;

    return 2;
}

#endif // SRC_HPP
