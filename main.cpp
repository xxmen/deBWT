#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <map>
#define DNA_LEN 5000000
#define KMER_LEN 20
#define KMER_SUM 4864385
using namespace std;
uint32_t get_c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    3};
char dna[DNA_LEN];
char bcode[DNA_LEN >> 2];
bool cmp(uint64_t a, uint64_t b)
{
    return (a << 2) < (b << 2);
}
const int kmer2_l = 2 * (KMER_LEN + 2);
int inputSeq(char *s, char *fName)
{

    FILE *fp;
    if((fp = fopen(fName, "r")) == NULL)
    {
        printf("文件不存在.\n");
        exit(1);
    }
    else
    {
        char info[100];
        fgets(info, 100, fp);
        uint32_t i = 0;
        while (~fscanf(fp, "%s", s + i * 70))
        {
            i++;
        }
        fclose(fp);
    }
    return (int)strlen(s);
}
uint64_t kc[KMER_SUM], K[KMER_SUM];
bool is_in[KMER_SUM], is_out[KMER_SUM];
uint64_t occ[KMER_SUM], ind[KMER_SUM];
size_t kid = 0;

int binSearch(uint64_t kmer)
{
    size_t pos;
    pos = lower_bound(K, K + kid, kmer) - K;
    return (int) (K[pos] == kmer ? pos : -1);
}

int main(int argc, const char * argv[])
{
    char fDNAname[] = "/Users/os/Desktop/deBWT/E.coli.fa";
    int dna_len = inputSeq(dna, fDNAname);
    FILE *fKmer;
    fKmer = fopen("/Users/os/Desktop/mer_counts_dumps.fa", "r");
    uint32_t cnt = 0;
    char kmer[22];
    size_t id = 0;
    while (~fscanf(fKmer, ">%d\n", &cnt))
    {
        uint64_t kmer_cnt = 0;
        fscanf(fKmer, "%s\n", kmer);
        for (int i = 0; i < KMER_LEN + 2; i++)
        {
            kmer_cnt = (kmer_cnt << 2) | get_c[kmer[i]];
        }
        uint64_t kcnt = cnt;
        kmer_cnt <<= (64 - 2 * (KMER_LEN + 2));
        kmer_cnt |= kcnt;
        kc[id++] = kmer_cnt;
    }
    sort(kc, kc + KMER_SUM, cmp);

    printf("%d %d\n", id, KMER_SUM);

    uint64_t mask_k = 1L, mask_l = -1L, mask_r = 1L, mask_n = 1L;
    size_t cnt_len = 64 - 2 * (KMER_LEN + 2);
    mask_k = ((mask_k << 62) - 1) >> (cnt_len + 2) << (cnt_len + 2);

    mask_l = (mask_l >> 62) << 62;

    mask_r = ((mask_r << (64 - 2 * (KMER_LEN + 2) + 2)) - 1) >> cnt_len << cnt_len;

    mask_n =  (mask_n << cnt_len) - 1;

    uint64_t x = mask_k;

    //cout<<bitset<64>(x)<<endl;

    size_t theindex = 0;

    uint32_t cnt_in = 0, cnt_out = 0;

    for (size_t i = 0, j, k; i < KMER_SUM; )
    {
        uint64_t now_kmer = kc[i] & mask_k;
        uint64_t pre = kc[i] & mask_l;
        uint64_t next = kc[i] & mask_r;

        bool multi_in = false;
        bool multi_out = false;
        for (j = i + 1; j < KMER_SUM; j++)
        {
            if ((kc[j] & mask_k) == now_kmer)
            {
                if (pre != (kc[j] & mask_l)) multi_in = true;
                if (next != (kc[j] & mask_r)) multi_out = true;
            }
            else
                break;
        }

        if (multi_in || multi_out)
        {
            if (multi_out)
            {
                is_out[kid] = true;
                cnt_out ++;
            }
            if (multi_in)
            {
                is_in[kid] = true;
                cnt_in ++;
                for (k = i; k < j; ++k)
                {
                    occ[kid] += (kc[k] & mask_n);
                }
                ind[kid] =  theindex;
                theindex += occ[kid];
            }
            K[kid] = (kc[i] & mask_k) >> (cnt_len + 2);
            ++kid;
        }
        i = j;
    }

    //for (int i = kid - 1; i >  kid - 20; i--)
     //   printf("%llu\n", K[i]);

    size_t BCN_size = (theindex);
    uint64_t *BCN;
    BCN = new uint64_t[BCN_size];

    uint64_t kmer_value = 0;
    uint64_t MASK = (1L << (2 * KMER_LEN)) - 1;
    size_t bcode_len = 1;
    uint64_t t_index = 0, tk;
    cout<<bitset<64>(MASK)<<endl;
    for (int i = 0, pos; i < dna_len; i++)
    {

        kmer_value = (kmer_value << 2) | get_c[dna[i]];
        kmer_value &= MASK;

        if (i >= KMER_LEN - 1)
        {
            //printf("%llu\n", kmer_value);
            pos = binSearch(kmer_value);
            //printf("%d\n", pos);
            if (pos != -1)
            {
                if (is_in[pos])  //multip in
                {
                    t_index = ind[pos];
                    tk = t_index;
                    while (tk < BCN_size && BCN[tk])
                        tk++;
                    if (i >= KMER_LEN)
                        BCN[tk] = (bcode_len << 3) | get_c[dna[i - KMER_LEN]];
                    else
                        BCN[tk] = (bcode_len << 3) | 4; //"$"
                }
                if (is_out[pos])
                {
                    bcode[bcode_len++] = dna[i + 1];
                }
            }

        }
    }

    bcode[bcode_len] = '\0';


    for (size_t i = 0, j, tmp_index=0, l_index=0, index; i<k2len;)
    {
        //find the io_info and check if mulip in
        j = i+1;
        uint64_t now_kmer = K2[i]&mask_c;
        while (l_index < kmer+1 && (last_string[l_index]&mask_c) <= now_kmer)
        {
            BWT[bwt_index++] = ((last_string[l_index++]&mask_l)>>62);
            // BWT[bwt_index-1] = 4;
        }
        if (now_kmer >= begin_kmer && begin_kmer >= last_kmer)
        {
            if (now_kmer == begin_kmer)
            {
                is_in = true;
            }
            else
            {
                BWT[bwt_index++] = 4;
            }
            //else if the begin_kmer no appear in K2, no need to create the io_info of begin_kmer
        }
        last_kmer = now_kmer;
        while (j<k2len && ((K2[j]&mask_c) == now_kmer))
        {
            if ((K2[j]&mask_l) != (K2[i]&mask_l)) is_in = true;
            if ((K2[j]&mask_r) !=( K2[i]&mask_r)) is_out = true;
            ++j;
        }
        if (is_in || is_out)
            tmp_mask = io_info[tmp_index++];
        if (is_in)
        {
            BCN_index = (tmp_mask&mask_index);
            BCN_len = ((tmp_mask&mask_in) >> 32);

            sort(BCN+BCN_index, BCN+BCN_index+BCN_len, bwt_cmp);

            for (size_t k=BCN_index; k<BCN_index+BCN_len; ++k)
            {
                BWT[bwt_index++] = BCN[k]&mask_code;
            }
        }
        else
        {
            tmp_num=0;
            for (size_t k = i; k<j; ++k)
            {
                tmp_num += (K2[k]&mask_n);
            }
            tmp_code = ((K2[i]&mask_l)>>62);
            while (tmp_num--)
            {
                BWT[bwt_index++] = tmp_code;
            }
        }
        i = j;
    }

    printf("%d\n", bcode_len);
    printf("%d   %d\n", cnt_in, cnt_out);
    printf("Time used = %.2f s\n",  (double)clock() / CLOCKS_PER_SEC);
    return 0;
}
