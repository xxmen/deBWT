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
uint8_t BWT[DNA_LEN + 10];
size_t bcode_len;
uint64_t MultiIn_INFO[10000];
bool cmp(uint64_t a, uint64_t b)
{
    return (a << 2) < (b << 2);
}
bool bcode_cmp(const uint64_t a, const uint64_t b)
{
    return (bcode[a >> 3] < bcode[b >> 3]);
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
uint64_t occ[KMER_SUM] = {0}, ind[KMER_SUM];
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
    fKmer = fopen("/Users/os/Desktop/deBWT/mer_counts_dumps.fa", "r");
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
    printf("Sort start = %.2f s\n",  (double)clock() / CLOCKS_PER_SEC);
    sort(kc, kc + KMER_SUM, cmp);

    printf("Sort finish = %.2f s\n",  (double)clock() / CLOCKS_PER_SEC);

    printf("%d %d\n", id, KMER_SUM);

    uint64_t mask_k = 1L, mask_l = -1L, mask_r = 1L, mask_n = 1L;
    size_t cnt_len = 64 - 2 * (KMER_LEN + 2);
    mask_k = ((mask_k << 62) - 1) >> (cnt_len + 2) << (cnt_len + 2);

    mask_l = (mask_l >> 62) << 62;

    mask_r = ((mask_r << (64 - 2 * (KMER_LEN + 2) + 2)) - 1) >> cnt_len << cnt_len;

    mask_n =  (mask_n << cnt_len) - 1;

    uint64_t x = mask_k;

    //cout<<bitset<64>(x)<<endl;

    uint32_t theindex = 0;

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
            for (k = i; k < j; ++k)
            {
                occ[kid] += (kc[k] & mask_n);
            }
            if (multi_out)
            {
                is_out[kid] = true;
                cnt_out ++;
            }
            if (multi_in)
            {
                is_in[kid] = true;
                cnt_in ++;
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

    uint32_t MultiIn_INFO_size = theindex;

    uint64_t kmer_value = 0;
    uint64_t MASK = (1L << (2 * KMER_LEN)) - 1;
    bcode_len = 1;
    uint64_t t_index = 0, tk;

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
                    while (tk < MultiIn_INFO_size && MultiIn_INFO[tk])
                        tk++;
                    if (i >= KMER_LEN)
                        MultiIn_INFO[tk] = (bcode_len << 3) | get_c[dna[i - KMER_LEN]];
                    else
                        MultiIn_INFO[tk] = (bcode_len << 3) | 4; //"$"
                }
                if (is_out[pos])
                {
                    bcode[bcode_len++] = dna[i + 1];
                }
            }

        }
    }

    bcode[bcode_len] = '\0';


    size_t tmp_index = 0;
    size_t bwt_index = 0;

    uint64_t mask_char = (1L << 3) - 1;
    uint64_t tmp = 0;

    for (size_t i = 0; i < KMER_LEN; ++i)
        tmp = (tmp << 2) | get_c[dna[i]];
    uint64_t begin_kmer = tmp << (64 - 2 * KMER_LEN) >> 2;

    uint64_t *last_string;
    last_string = new uint64_t[KMER_LEN + 1];
    tmp = 0;
    for (size_t i = KMER_LEN + 1; i >= 1; --i)
    {
        for (size_t j = dna_len - i; j < dna_len; ++j)
            tmp = (tmp << 2) | get_c[dna[j]];
        last_string[KMER_LEN + 1 - i] = tmp << (64 - (i) * 2);
    }
    sort(last_string, last_string + KMER_LEN + 1, cmp);

    uint64_t last_kmer = 0;


    int Len2 = 0, Len1 = 0, Len3 = 0;
    for (size_t i = 0, j, l_index = 0, index; i < KMER_SUM; )
    {
        uint64_t now_kmer = kc[i] & mask_k;

        while (l_index < KMER_LEN + 1 && (last_string[l_index] & mask_k) <= now_kmer)
            BWT[bwt_index++] = ((last_string[l_index++] & mask_l) >> 62);

        if (now_kmer >= begin_kmer && begin_kmer >= last_kmer)
        {
            BWT[bwt_index++] = 4;
        }

        last_kmer = now_kmer;

        if (now_kmer != (kc[i + 1] & mask_k))    //Case 1
        {
            int tmp_cnt = kc[i] & mask_n;
            Len1 += tmp_cnt;
            uint8_t tmp_code = ((kc[i] & mask_l) >> 62);
            while (tmp_cnt--)
                BWT[bwt_index++] = tmp_code;
            i++;
            continue;
        }
        if (is_in[tmp_index])                           //Case 3
        {
            Len3 += occ[tmp_index];
            sort(MultiIn_INFO + ind[tmp_index], MultiIn_INFO + ind[tmp_index] + occ[tmp_index], bcode_cmp);
            for (size_t k = ind[tmp_index]; k < ind[tmp_index] + occ[tmp_index]; ++k)
            {
                BWT[bwt_index++] = MultiIn_INFO[k] & mask_char;
            }
        }
        else                                           //Case 2
        {
            size_t tmp_cnt = occ[tmp_index];
            Len2 += tmp_cnt;
            uint8_t tmp_code = ((kc[i] & mask_l) >> 62);
            while (tmp_cnt--)
            {
                BWT[bwt_index++] = tmp_code;
            }
        }
        tmp_index++;

        for (j = i + 1; j < KMER_SUM && (kc[j] & mask_k) == now_kmer; )
            j++;
        i = j;
    }
    printf("len1 = %d  len2 = %d len3 = %d\n", Len1, Len2, Len3);
    printf("bcode_len = %d\n", bcode_len);
    printf("%d   %d\n", cnt_in, cnt_out);

    cout<<bwt_index<<"  "<<dna_len<<endl;
    FILE *fout = fopen("/Users/os/Desktop/deBWT/res.bwt","w");
    char cc[6]={'A', 'C', 'G', 'T', '$', 'X'};
    for (size_t i = 0; i <= dna_len; ++i)
    {
        fprintf(fout, "%c", cc[BWT[i]]);
        if ((i+1) % 80 == 0)
            fprintf(fout, "\n");
    }

    printf("Time used = %.2f s\n",  (double)clock() / CLOCKS_PER_SEC);
    return 0;
}
