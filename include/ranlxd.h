#ifndef INCLUDED_RANLUX
#define INCLUDED_RANLUX
#include <complex>

// This is Luscher's ranlux random generator which I turned into a C++ version
// I've also added code for computing Gaussian distributions using the Box-Mueller method
// A global random generator is instantiated at the end of the file, called "rng"

//---- Ranlxd random generator declarations
typedef struct
{
   int c1,c2,c3,c4;
} vec_t;

typedef struct
{
   vec_t c1,c2;
} dble_vec_t;

class Ranlxd
{
  private:
//    int init;
    int pr;
    int prm;
    int ir;
    int jr;
    int is;
    int is_old;
    int next[96];
    float one_bit;
    vec_t carry;
    union { dble_vec_t vec[12]; int num[96]; } x;

    void Error(int);
    void Update(void);
    void DefineConstants(void);
    void LuxLevel(int);

    int lux_level;

    Ranlxd& operator=(const Ranlxd&);
    Ranlxd(const Ranlxd&);

  public:
    explicit Ranlxd(int);
    ~Ranlxd();
    void Init(long);
    double UniformDeviate();
    double Normal(double, double);
    std::complex<double> Normal(double);
    int StateSize();
    void GetState(long*);
    void SetState(long*);


};

extern Ranlxd rng; // global generator initialized with lux-level = 2

#endif
