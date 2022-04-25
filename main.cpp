/* Copyright (C) 2019-2021 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

// This is a sample program for education purposes only.
// It attempts to show the various basic mathematical
// operations that can be performed on both ciphertexts
// and plaintexts.

#include <iostream>
#include "myutils.h"
#include <helib/helib.h>
#include <thread>

//#define DBG_INFO

// NOTE: reference: 'tests/GTestThinBootstrapping.cpp'

constexpr static long mValues[][14] = {
        // clang-format off
        //{  p,phi(m),     m,  d,  m1,  m2, m3,    g1,    g2,   g3,ord1,ord2,ord3, c_m}
        {  2,    48,   105, 12,   3,  35,  0,    71,    76,    0,   2,   2,   0, 100},
        {  2,   600,  1023, 10,  11,  93,  0,   838,   584,    0,  10,   6,   0, 100}, // m=(3)*11*{31} m/phim(m)=1.7    C=24  D=2 E=1
        {  2,  1200,  1705, 20,  11, 155,  0,   156,   936,    0,  10,   6,   0, 100}, // m=(5)*11*{31} m/phim(m)=1.42   C=34  D=2 E=2
        {  2,  1728,  4095, 12,   7,   5,117,  2341,  3277, 3641,   6,   4,   6, 100}, // m=(3^2)*5*7*{13} m/phim(m)=2.36 C=26 D=3 E=2
        {  2,  2304,  4641, 24,   7,   3,221,  3979,  3095, 3760,   6,   2,  -8, 300}, // m=3*7*(13)*{17} :-( m/phim(m)=2.01 C=45 D=4 E=3
        {  2,  4096,  4369, 16,  17, 257,  0,   258,  4115,    0,  16, -16,   0, 100}, // m=17*(257) :-( m/phim(m)=1.06 C=61 D=3 E=4
        {  2, 12800, 17425, 40,  41, 425,  0,  5951,  8078,    0,  40,  -8,   0, 100}, // m=(5^2)*{17}*41 m/phim(m)=1.36 C=93  D=3 E=3
        {  2, 15004, 15709, 22,  23, 683,  0,  4099, 13663,    0,  22,  31,   0, 100}, // m=23*(683) m/phim(m)=1.04      C=73  D=2 E=1
        {  2, 16384, 21845, 16,  17,   5,257,  8996, 17477,21591,  16,   4, -16, 200}, // m=5*17*(257) :-( m/phim(m)=1.33 C=65 D=4 E=4
        {  2, 18000, 18631, 25,  31, 601,  0, 15627,  1334,    0,  30,  24,   0, 100}, // m=31*(601) m/phim(m)=1.03      C=77  D=2 E=0
        {  2, 18816, 24295, 28,  43, 565,  0, 16386, 16427,    0,  42,  16,   0, 100}, // m=(5)*43*{113} m/phim(m)=1.29  C=84  D=2 E=2
        {  2, 21168, 27305, 28,  43, 635,  0, 10796, 26059,    0,  42,  18,   0, 100}, // m=(5)*43*{127} m/phim(m)=1.28  C=86  D=2 E=2
        {  2, 23040, 28679, 24,  17,   7,241, 15184,  4098,28204,  16,   6, -10, 200}, // m=7*17*(241) m/phim(m)=1.24    C=63  D=4 E=3
        {  2, 24000, 31775, 20,  41, 775,  0,  6976, 24806,    0,  40,  30,   0, 100}, // m=(5^2)*{31}*41 m/phim(m)=1.32 C=88  D=2 E=2
        {  2, 26400, 27311, 55,  31, 881,  0, 21145,  1830,    0,  30,  16,   0, 100}, // m=31*(881) m/phim(m)=1.03      C=99  D=2 E=0
        {  2, 27000, 32767, 15,  31,   7,151, 11628, 28087,25824,  30,   6, -10, 200},
        {  2, 31104, 35113, 36,  37, 949,  0, 16134,  8548,    0,  36,  24,   0, 200}, // m=(13)*37*{73} m/phim(m)=1.12  C=94  D=2 E=2
        {  2, 34848, 45655, 44,  23,1985,  0, 33746, 27831,    0,  22,  36,   0, 100}, // m=(5)*23*{397} m/phim(m)=1.31  C=100 D=2 E=2
        {  2, 42336, 42799, 21, 127, 337,  0, 25276, 40133,    0, 126,  16,   0, 200}, // m=127*(337) m/phim(m)=1.01     C=161 D=2 E=0
        {  2, 45360, 46063, 45,  73, 631,  0, 35337, 20222,    0,  72,  14,   0, 100}, // m=73*(631) m/phim(m)=1.01      C=129 D=2 E=0
        {  2, 46080, 53261, 24,  17,  13,241, 43863, 28680,15913,  16,  12, -10, 100}, // m=13*17*(241) m/phim(m)=1.15   C=69  D=4 E=3
        {  2, 49500, 49981, 30, 151, 331,  0,  6952, 28540,    0, 150,  11,   0, 100}, // m=151*(331) m/phim(m)=1        C=189 D=2 E=1
        {  2, 54000, 55831, 25,  31,1801,  0, 19812, 50593,    0,  30,  72,   0, 100}, // m=31*(1801) m/phim(m)=1.03     C=125 D=2 E=0
        {  2, 60016, 60787, 22,  89, 683,  0,  2050, 58741,    0,  88,  31,   0, 200}, // m=89*(683) m/phim(m)=1.01      C=139 D=2 E=1

        {  7,    36,    57,  3,   3,  19,  0,    20,    40,    0,   2,  -6,   0, 100}, // m=3*(19) :-( m/phim(m)=1.58 C=14 D=3 E=0

        { 17,    48,   105, 12,   3,  35,  0,    71,    76,    0,   2,   2,   0, 100}, // m=3*(5)*{7} m/phim(m)=2.18 C=14 D=2 E=2
        { 17,   576,  1365, 12,   7,   3, 65,   976,   911,  463,   6,   2,   4, 100}, // m=3*(5)*7*{13} m/phim(m)=2.36  C=22  D=3
        { 17, 18000, 21917, 30, 101, 217,  0,  5860,  5455,    0, 100,   6,   0, 100}, // m=(7)*{31}*101 m/phim(m)=1.21  C=134 D=2
        { 17, 30000, 34441, 30, 101, 341,  0,  2729, 31715,    0, 100,  10,   0, 100}, // m=(11)*{31}*101 m/phim(m)=1.14 C=138 D=2
        { 17, 40000, 45551, 40, 101, 451,  0, 19394,  7677,    0, 100,  10,   0, 200}, // m=(11)*{41}*101 m/phim(m)=1.13 C=148 D=2
        { 17, 46656, 52429, 36, 109, 481,  0, 46658,  5778,    0, 108,  12,   0, 100}, // m=(13)*{37}*109 m/phim(m)=1.12 C=154 D=2
        { 17, 54208, 59363, 44,  23,2581,  0, 25811,  5199,    0,  22,  56,   0, 100}, // m=23*(29)*{89} m/phim(m)=1.09  C=120 D=2
        { 17, 70000, 78881, 10, 101, 781,  0, 67167, 58581,    0, 100,  70,   0, 100}, // m=(11)*{71}*101 m/phim(m)=1.12 C=178 D=2

        {127,   576,  1365, 12,  7,    3, 65,   976,   911,  463,   6,   2,   4, 100}, // m=3*(5)*7*{13} m/phim(m)=2.36   C=22  D=3
        {127,  1200,  1925, 20, 11,  175,  0,  1751,   199,    0,  10,   6,   0, 100}, //  m=(5^2)*{7}*11 m/phim(m)=1.6   C=34 D=2
        {127,  2160,  2821, 30, 13,  217,  0,   652,   222,    0,  12,   6,   0, 100}, // m=(7)*13*{31} m/phim(m)=1.3     C=46 D=2
        {127, 18816, 24295, 28, 43,  565,  0, 16386, 16427,    0,  42,  16,   0, 100}, // m=(5)*43*{113} m/phim(m)=1.29   C=84  D=2
        {127, 26112, 30277, 24, 17, 1781,  0, 14249, 10694,    0,  16,  68,   0, 100}, // m=(13)*17*{137} m/phim(m)=1.15  C=106 D=2
        {127, 31752, 32551, 14, 43,  757,  0,  7571, 28768,    0,  42,  54,   0, 100}, // m=43*(757) :-( m/phim(m)=1.02   C=161 D=3
        {127, 46656, 51319, 36, 37, 1387,  0, 48546, 24976,    0,  36, -36,   0, 200}, // m=(19)*37*{73}:-( m/phim(m)=1.09 C=141 D=3
        {127, 49392, 61103, 28, 43, 1421,  0,  1422, 14234,    0,  42,  42,   0, 200}, // m=(7^2)*{29}*43 m/phim(m)=1.23  C=110 D=2
        {127, 54400, 61787, 40, 41, 1507,  0, 30141, 46782,    0,  40,  34,   0, 100}, // m=(11)*41*{137} m/phim(m)=1.13  C=112 D=2
        {127, 72000, 77531, 30, 61, 1271,  0,  7627, 34344,    0,  60,  40,   0, 100}  // m=(31)*{41}*61 m/phim(m)=1.07   C=128 D=2
        // clang-format on
};

constexpr static long num_mValues = sizeof(mValues) / (14 * sizeof(long));

static std::size_t getIdx(const long p, const long N)
{
    for (long i = 0; i < (long)num_mValues; i++)
        if (mValues[i][0] == p && mValues[i][1] >= N) {
            return i;
        }
    throw std::invalid_argument(
            "Could not find mValues row corresponding to p=" + std::to_string(p) +
            ", N=" + std::to_string(N));
}

static long validatem(const long m, const long p)
{
    if (NTL::GCD(p, m) != 1)
        throw std::invalid_argument("m and p must be coprime");
    return m;
}

static NTL::Vec<long> calculateMvec(const std::size_t idx)
{
    NTL::Vec<long> mvec;
    append(mvec, mValues[idx][4]);
    if (mValues[idx][5] > 1)
        append(mvec, mValues[idx][5]);
    if (mValues[idx][6] > 1)
        append(mvec, mValues[idx][6]);
    return mvec;
}

static std::vector<long> calculateGens(const std::size_t idx)
{
    std::vector<long> gens;
    gens.push_back(mValues[idx][7]);
    if (mValues[idx][8] > 1)
        gens.push_back(mValues[idx][8]);
    if (mValues[idx][9] > 1)
        gens.push_back(mValues[idx][9]);
    return gens;
}

static std::vector<long> calculateOrds(const std::size_t idx)
{
    std::vector<long> ords;
    ords.push_back(mValues[idx][10]);
    if (std::abs(mValues[idx][11]) > 1)
        ords.push_back(mValues[idx][11]);
    if (std::abs(mValues[idx][12]) > 1)
        ords.push_back(mValues[idx][12]);
    return ords;
}

struct noise_stat{
    double cap, noise, mod;
    double mean, var, mag;
};

noise_stat decrypt_stat(const helib::Ctxt& ctxt, const helib::Ptxt<helib::BGV>& expected_ptxt,
                        const helib::SecKey& secKey, int curlevel, NTL::ZZX* dst=nullptr){
    auto& context = ctxt.getContext();

    NTL::ZZX ptxt_poly = expected_ptxt.getPolyRepr();
    // NOTE: helib::PolyRed will keep the sign of coefficients when q=2
//    helib::PolyRed(ptxt_poly, ctxt.getPtxtSpace());
    auto q = ctxt.getPtxtSpace();
    ZZX_mod(ptxt_poly, q);


    NTL::ZZX _, poly;
    secKey.Decrypt(_, ctxt, poly);
    helib::PolyRed(poly, context.productOfPrimes(ctxt.getPrimeSet()));
    poly.normalize();

    if(dst){
        *dst = poly;
        ZZX_mod(*dst, q);
        dst->normalize();
    }

    double mean, var, mag;
    double cap = ctxt.capacity(),
        noise = std::log2(to_double(ctxt.getNoiseBound())),
        mod = std::log2(NTL::to_double(context.productOfPrimes(ctxt.getPrimeSet())));
    mean_var_mag(poly - ptxt_poly, mean, var, mag, context.getPhiM());
#ifdef DBG_INFO
    std::cout << "level = " << curlevel <<
              ", capacity = " << cap <<
              ", noise = " << noise <<
              // NOTE: context.logOfPrimeSet compute the NATURAL logarithm
              ", modulus = " << mod <<
              ", noise mean = " << mean <<
              ", noise var = " << var <<
              ", noise mag = " << mag <<
              '\n';
#endif
    return {cap, noise, mod, mean, var, mag};
}

void test(const helib::Context* c){}

// worker thread
void thread_func(int param_idx, int samples_count, const std::string& output_file,
                 int total_levels, int boot_levels, int square_after_boot,
                 const helib::Context& context, const helib::PubKey& pubKey,
                 const helib::SecKey& secKey){
    auto f = std::fopen(output_file.c_str(), "w");
    fprintf(f, "%d %d %d %d\n", param_idx, total_levels, boot_levels, square_after_boot);

    for(int sample = 0; sample < samples_count; sample++) {
//        std::cout << "sample " << sample << '\n';

        helib::Ptxt<helib::BGV> ptxt(context);
        for (int i = 0; i < ptxt.size(); ++i) {
            ptxt[i] = 1;
        }

        helib::Ctxt ctxt(pubKey);
        pubKey.Encrypt(ctxt, ptxt);

        std::vector<noise_stat> noise_stat_vec;

        int levels = 0;
        while (ctxt.capacity() > 0 && levels < total_levels) {

            noise_stat_vec.push_back(decrypt_stat(ctxt, ptxt, secKey, levels));

            ctxt.square();
            ptxt.square(); // NOTE: ptxt stores only the slots, not the polynomial (which is a natural idea for efficiency)
            levels++;
        }

        noise_stat_vec.push_back(decrypt_stat(ctxt, ptxt, secKey, levels));

        for(int blvl = 0; blvl < boot_levels; blvl++) {
#ifdef DBG_INFO
            std::cout << "start bootstrapping No." << blvl << "\n";
#endif
            pubKey.thinReCrypt(ctxt);

            // Decrypt the modified ciphertext into a new plaintext
            NTL::ZZX plaintext_poly_result;
            noise_stat_vec.push_back(decrypt_stat(ctxt, ptxt, secKey, -1, &plaintext_poly_result));


            auto real_poly = ptxt.getPolyRepr();
            ZZX_mod(real_poly, ctxt.getPtxtSpace());
            if (plaintext_poly_result != real_poly) {
                std::cout << deg(plaintext_poly_result) << " xxx " << deg(real_poly) << '\n';
                for (int i = 0, end = std::min(10, (int) std::min(deg(plaintext_poly_result), deg(real_poly)));
                     i < end; i++) {
                    std::cout << plaintext_poly_result[i] << " # " << real_poly[i] << '\n';
                }
                std::cout << "ERROR: decryption result doesn't match expected value\n";
            } else {
#ifdef DBG_INFO
                std::cout << "Decryption ok\n";
#endif
            }

            if(blvl == boot_levels - 1) // jump out of loop right after the last bootstrapping
                break;

            for(int slvl = 0; slvl < square_after_boot; slvl++) {
                // one more square
                ctxt.square();
                ptxt.square();
                noise_stat_vec.push_back(decrypt_stat(ctxt, ptxt, secKey, -2-slvl));
            }
        }

        // NOTE: noise_stat_vec has the stat of level 0, 1, ... total_levels, total_levels+1(after bootstrapping)
        fprintf(f, "SAMPLE %d\n", sample);
        for(auto stat : noise_stat_vec){
            fprintf(f, "%f, %f, %f, %f, %f, %f\n",
                    stat.cap, stat.noise, stat.mod,
                    stat.mean, stat.var, stat.mag);
        }
        std::fflush(f);
    }
    std::fclose(f);
}

int main(int argc, char* argv[]) {
    /*  Example of BGV scheme  */

//    int test_idx;
//    int test_r;
//    int test_bits;
//    while(true){
//        std::cout << "input idx:";
//        std::cin >> test_idx;
//        std::cout << "input r:";
//        std::cin >> test_r;
//        std::cout << "input bits:";
//        std::cin >> test_bits;
//        if(test_idx < 0 || test_idx >= num_mValues || test_r <= 0 || test_bits <= 0)
//            return 0;
//
//        auto p = mValues[test_idx][0];
//        auto m = mValues[test_idx][2];
//        auto mvec = calculateMvec(test_idx);
//        auto gens = calculateGens(test_idx);
//        auto ords = calculateOrds(test_idx);
//
//        helib::Context context = helib::ContextBuilder<helib::BGV>()
//                .m(m)
//                .p(p)
//                .bits(test_bits)
//                .r(test_r)
//                .c(3)
//                .mvec(mvec)
//                .ords(ords)
//                .gens(gens)
//                .bootstrappable()
//                .thinboot()
//                .buildCache(true)
//                .build();
//        std::cout << "security level: " << context.securityLevel() << '\n';
//        std::cout << "prime size: " << context.logOfPrime(context.getCtxtPrimes(1).first()) << '\n';
//        std::cout << "q size: " << context.bitSizeOfQ() << '\n';
//    }

    /**
     * idx  bits    valid
     * 7    180     x
     * 8    235     x
     * 9    275     1
     * 10   285     1
     * 11   340     v(1)
     * 12   359     v(2)
     */

    // Plaintext prime modulus
    unsigned long p = 2;
    // Hensel lifting (default = 1) NOTE: plaintext space is p^r
    unsigned long r = 1;
    // Number of columns of Key-Switching matrix (default = 2 or 3)
    unsigned long c = 3;
    // Number of bits of the modulus chain
    unsigned long bits = 359;
    // NOTE: least phi(m)
    unsigned long N = 512;


    bool force_idx = true;
    int param_idx = 12;
    // NOTE: parameters from hom aes
    if(!force_idx)
        param_idx = getIdx(p, N);
    p = mValues[param_idx][0];

    std::cout << "using parameter set No." << param_idx << '\n';

    // Cyclotomic polynomial - defines phi(m)
    unsigned long m = mValues[param_idx][2];
    auto mvec = calculateMvec(param_idx);
    auto gens = calculateGens(param_idx);
    auto ords = calculateOrds(param_idx);

    // NOTE: use the parameters from homomorphic aes


    std::cout << "\n*********************************************************";
    std::cout << "\n*         Basic Mathematical Operations Example         *";
    std::cout << "\n*         =====================================         *";
    std::cout << "\n*                                                       *";
    std::cout << "\n* This is a sample program for education purposes only. *";
    std::cout << "\n* It attempts to show the various basic mathematical    *";
    std::cout << "\n* operations that can be performed on both ciphertexts  *";
    std::cout << "\n* and plaintexts.                                       *";
    std::cout << "\n*                                                       *";
    std::cout << "\n*********************************************************";
    std::cout << std::endl;

    std::cout << "Initialising context object..." << std::endl;
    // Initialize context
    // This object will hold information about the algebra created from the
    // previously set parameters
    helib::Context context = helib::ContextBuilder<helib::BGV>()
            .m(m)
            .p(p)
            .bits(bits)
            .r(r)
            .c(c)
            .mvec(mvec)
            .ords(ords)
            .gens(gens)
            .bootstrappable()
            .thinboot()
            .buildCache(true)
            .build();

    // Print the context
    context.printout();
    std::cout << std::endl;
    // NOTE: print prime bit lengths
    std::cout << "prime lengths: ";
    for(auto i : context.getCtxtPrimes()){
        std::cout << i << "->" << context.logOfPrime(i) << " ";
    }
    std::cout << '\n';

    // Print the security level
    std::cout << "Security: " << context.securityLevel() << std::endl;

    // Secret key management
    std::cout << "Creating secret key..." << std::endl;
    // Create a secret key associated with the context
    // NOTE: enableBootStrapping has no effect because in it is already called in contextBuilder
//    context.enableBootStrapping(helib::convert<NTL::Vec<long>>(mvec));
    helib::SecKey secret_key(context);
    // Generate the secret key
    secret_key.GenSecKey();
    // Compute key-switching matrices that we need
    helib::addSome1DMatrices(secret_key);
    helib::addFrbMatrices(secret_key);
    secret_key.genRecryptData();
    std::cout << "Generating key-switching matrices..." << std::endl;

    // Public key management
    // Set the secret key (upcast: SecKey is a subclass of PubKey)
    helib::PubKey public_key = secret_key; // NOTE: when using copy constructor, the public key is a real 'public' key


    // Get the EncryptedArray of the context
    const helib::EncryptedArray &ea = context.getEA();

    // Get the number of slot (phi(m))
    long nslots = ea.size();
    std::cout << "Number of slots: " << nslots << std::endl;

    double level_estimate = context.logOfProduct(context.getCtxtPrimes()) /
            context.logOfPrime(context.getCtxtPrimes(1).first());
    std::cout << "Number of level(estimated): " << level_estimate << '\n';

    /**
     * ############################
     * Get lots of samples
     * ############################
     */
    int samples_count = 100000; // something that will never be reached...
    std::string output_file = argc > 1 ? argv[1] : "/dev/null";
    int total_levels = 19;
    int boot_levels = 2;
    int square_after_boot = 2;

    int workers = argc > 2 ? std::stoi(argv[2]) : 1; // default to 1 thread
    std::vector<std::thread> threads;
    threads.reserve(workers);
    for(int i = 0; i < workers; i++){
        threads.emplace_back(thread_func,
                                      param_idx,
                                      samples_count,
                                      output_file + "_" + std::to_string(i),
                                      total_levels,
                                      boot_levels,
                                      square_after_boot,
                                      std::cref(context),
                                      std::cref(public_key),
                                      std::cref(secret_key));
    }

    for(auto& t : threads)
        t.join();

    return 0;
}