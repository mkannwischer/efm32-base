/// @file rainbow.c
/// @brief The standard implementations for functions in rainbow.h
///
#include "rainbow_config.h"

#include "rainbow_keypair.h"

#include "rainbow.h"

#include "blas.h"

#include "rainbow_blas.h"
#include "rainbow_asm.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "utils_prng.h"
#include "utils_hash.h"
#include "utils_malloc.h"


#define MAX_ATTEMPT_FRMAT  128


/////////////////////////////


#define _MAX_O  ((_O1>_O2)?_O1:_O2)
#define _MAX_O_BYTE  ((_O1_BYTE>_O2_BYTE)?_O1_BYTE:_O2_BYTE)



int rainbow_sign( uint8_t * signature , const sk_t * sk , const uint8_t * _digest )
{
    // allocate temporary storage.
    uint8_t mat_l1[_O1*_O1_BYTE];
    uint8_t mat_l2[_O2*_O2_BYTE];
    uint8_t mat_buffer[2*_MAX_O*_MAX_O_BYTE];

    // setup PRNG
    prng_t prng_sign;
    uint8_t prng_preseed[LEN_SKSEED+_HASH_LEN];
    memcpy( prng_preseed , sk->sk_seed , LEN_SKSEED );
    memcpy( prng_preseed + LEN_SKSEED , _digest , _HASH_LEN );                        // prng_preseed = sk_seed || digest
    uint8_t prng_seed[_HASH_LEN];
    hash_msg( prng_seed , _HASH_LEN , prng_preseed , _HASH_LEN+LEN_SKSEED );
    prng_set( &prng_sign , prng_seed , _HASH_LEN );                                   // seed = H( sk_seed || digest )
    for(unsigned i=0;i<LEN_SKSEED+_HASH_LEN;i++) prng_preseed[i] ^= prng_preseed[i];  // clean
    for(unsigned i=0;i<_HASH_LEN;i++) prng_seed[i] ^= prng_seed[i];                   // clean

    // roll vinegars.
    uint8_t vinegar[_V1_BYTE];
    unsigned n_attempt = 0;
    unsigned l1_succ = 0;
    while( !l1_succ ) {
        if( MAX_ATTEMPT_FRMAT <= n_attempt ) break;
        prng_gen( &prng_sign , vinegar , _V1_BYTE );                       // generating vinegars
        #if ( _O1*_O1_BYTE  == 512) && (_V1 == 32)
        //gfmat_prod( mat_l1 , sk->l1_F2 , _O1*_O1_BYTE , _V1 , vinegar );   // generating the linear equations for layer 1
        gf16mat_prod_512_32_asm( mat_l1 , sk->l1_F2, vinegar );   // generating the linear equations for layer 1
        #elif ( _O1*_O1_BYTE  == 512) && (_V1 == 36)
        gf16mat_prod_512_36_asm( mat_l1 , sk->l1_F2, vinegar );   // generating the linear equations for layer 1
        #else
            #error not implemented
        #endif
        l1_succ = gfmat_inv( mat_l1 , mat_l1 , _O1 , mat_buffer );         // check if the linear equation solvable
        n_attempt ++;
    }

    // Given the vinegars, pre-compute variables needed for layer 2
    uint8_t r_l1_F1[_O1_BYTE] = {0};
    uint8_t r_l2_F1[_O2_BYTE] = {0};
    batch_quad_trimat_eval( r_l1_F1, sk->l1_F1, vinegar, _V1, _O1_BYTE );
    batch_quad_trimat_eval( r_l2_F1, sk->l2_F1, vinegar, _V1, _O2_BYTE );
    uint8_t mat_l2_F3[_O2*_O2_BYTE];
    uint8_t mat_l2_F2[_O1*_O2_BYTE];

    #if (_O2*_O2_BYTE == 512) && (_V1 == 32)
    //gfmat_prod( mat_l2_F3 , sk->l2_F3 , _O2*_O2_BYTE , _V1 , vinegar );
    gf16mat_prod_512_32_asm( mat_l2_F3 , sk->l2_F3 , vinegar );
    #elif (_O2*_O2_BYTE == 512) && (_V1 == 36)
    gf16mat_prod_512_36_asm( mat_l2_F3 , sk->l2_F3 , vinegar );
    #else
        #error not implemented
    #endif
    #if (_O1*_O2_BYTE == 512) && (_V1 == 32)
    //gfmat_prod( mat_l2_F2 , sk->l2_F2 , _O1*_O2_BYTE , _V1 , vinegar );
    gf16mat_prod_512_32_asm( mat_l2_F2 , sk->l2_F2, vinegar );
    #elif (_O1*_O2_BYTE == 512) && (_V1 == 36)
    gf16mat_prod_512_36_asm( mat_l2_F2 , sk->l2_F2, vinegar );
    #else
        #error not implemented
    #endif

    // Some local variables.
    uint8_t _z[_PUB_M_BYTE];
    uint8_t y[_PUB_M_BYTE];
    uint8_t * x_v1 = vinegar;
    uint8_t x_o1[_O1_BYTE];
    uint8_t x_o2[_O1_BYTE];

    uint8_t digest_salt[_HASH_LEN + _SALT_BYTE];
    memcpy( digest_salt , _digest , _HASH_LEN );
    uint8_t * salt = digest_salt + _HASH_LEN;

    uint8_t temp_o[_MAX_O_BYTE + 32]  = {0};
    unsigned succ = 0;
    while( !succ ) {
        if( MAX_ATTEMPT_FRMAT <= n_attempt ) break;
        // The computation:  H(digest||salt)  -->   z   --S-->   y  --C-map-->   x   --T-->   w

        prng_gen( &prng_sign , salt , _SALT_BYTE );                        // roll the salt
        hash_msg( _z , _PUB_M_BYTE , digest_salt , _HASH_LEN+_SALT_BYTE ); // H(digest||salt)

        //  y = S^-1 * z
        memcpy(y, _z, _PUB_M_BYTE);                                     // identity part of S
        #if (_O1_BYTE == 16) && (_O2 == 32)
        //gfmat_prod(temp_o, sk->s1, _O1_BYTE, _O2, _z+_O1_BYTE);
        gf16mat_prod_16_32_asm(temp_o, sk->s1, _z+_O1_BYTE);
        #else
            #error not implemented
        #endif
        gf256v_add(y, temp_o, _O1_BYTE);

        // Central Map:
        // layer 1: calculate x_o1
        memcpy( temp_o , r_l1_F1 , _O1_BYTE );
        gf256v_add( temp_o , y , _O1_BYTE );
        #if (_O1_BYTE == 16) && (_O1 == 32)
        //gfmat_prod( x_o1 , mat_l1, _O1_BYTE , _O1 , temp_o );
        gf16mat_prod_16_32_asm( x_o1 , mat_l1, temp_o );
        #else
            #error not implemented
        #endif
        // layer 2: calculate x_o2
        gf256v_set_zero( temp_o , _O2_BYTE );
        #if (_O2_BYTE == 16) && (_O1 == 32)
        //gfmat_prod( temp_o , mat_l2_F2, _O2_BYTE , _O1 , x_o1 );            // F2
        gf16mat_prod_16_32_asm( temp_o , mat_l2_F2, x_o1 );            // F2
        #else
            #error not implemented
        #endif

        batch_quad_trimat_eval( mat_l2 , sk->l2_F5, x_o1 , _O1, _O2_BYTE ); // F5
        gf256v_add( temp_o , mat_l2 , _O2_BYTE );
        gf256v_add( temp_o , r_l2_F1 , _O2_BYTE );                      // F1
        gf256v_add( temp_o , y + _O1_BYTE , _O2_BYTE );

        // generate the linear equations of the 2nd layer
        #if (_O2*_O2_BYTE  == 512) && (_O1 == 32)
        //gfmat_prod( mat_l2 , sk->l2_F6 , _O2*_O2_BYTE , _O1 , x_o1 );   // F6
        gf16mat_prod_512_32_asm( mat_l2 , sk->l2_F6 , x_o1 );   // F6
        #else
            #error not implemented
        #endif
        gf256v_add( mat_l2 , mat_l2_F3 , _O2*_O2_BYTE);                 // F3
        succ = gfmat_inv( mat_l2 , mat_l2 , _O2 , mat_buffer );
        #if (_O2_BYTE == 16) && (_O2 == 32)
        //gfmat_prod( x_o2 , mat_l2 , _O2_BYTE , _O2 , temp_o );         // solve l2 eqs
        gf16mat_prod_16_32_asm( x_o2 , mat_l2 , temp_o );         // solve l2 eqs
        #else
            #error not implemented
        #endif


        n_attempt ++;
    };
    //  w = T^-1 * y
    uint8_t w[_PUB_N_BYTE];
    // identity part of T.
    memcpy( w , x_v1 , _V1_BYTE );
    memcpy( w + _V1_BYTE , x_o1 , _O1_BYTE );
    memcpy( w + _V2_BYTE , x_o2 , _O2_BYTE );
    // Computing the t1 part.
    #if (_V1_BYTE == 16) && (_O1 == 32)
    //gfmat_prod(y, sk->t1, _V1_BYTE , _O1 , x_o1 );
    gf16mat_prod_16_32_asm(y, sk->t1, x_o1 );
    #elif (_V1_BYTE == 18) && (_O1 == 32)
    gf16mat_prod_18_32_asm(y, sk->t1, x_o1 );
    #else
        #error not implemented
    #endif
    gf256v_add(w, y, _V1_BYTE );
    // Computing the t4 part.
    #if (_V1_BYTE == 16) && (_O2 == 32)
    //gfmat_prod(y, sk->t4, _V1_BYTE , _O2 , x_o2 );
    gf16mat_prod_16_32_asm(y, sk->t4, x_o2 );
    #elif (_V1_BYTE == 18) && (_O1 == 32)
    gf16mat_prod_18_32_asm(y, sk->t4, x_o2 );
    #else
        #error not implemented
    #endif
    gf256v_add(w, y, _V1_BYTE );
    // Computing the t3 part.
    #if (_O1_BYTE == 16) && (_O2 == 32)
    //gfmat_prod(y, sk->t3, _O1_BYTE , _O2 , x_o2 );
    gf16mat_prod_16_32_asm(y, sk->t3, x_o2 );
    #else
        #error not implemented
    #endif
    gf256v_add(w+_V1_BYTE, y, _O1_BYTE );

    memset( signature , 0 , _SIGNATURE_BYTE );  // set the output 0
    // clean
    memset( mat_l1 , 0 , _O1*_O1_BYTE );
    memset( mat_l2 , 0 , _O2*_O2_BYTE );
    memset( mat_buffer , 0 , 2*_MAX_O*_MAX_O_BYTE );
    memset( &prng_sign , 0 , sizeof(prng_t) );
    memset( vinegar , 0 , _V1_BYTE );
    memset( r_l1_F1 , 0 , _O1_BYTE );
    memset( r_l2_F1 , 0 , _O2_BYTE );
    memset( mat_l2_F3 , 0 , _O2*_O2_BYTE );
    memset( mat_l2_F2 , 0 , _O1*_O2_BYTE );
    memset( _z , 0 , _PUB_M_BYTE );
    memset( y , 0 , _PUB_M_BYTE );
    memset( x_o1 , 0 , _O1_BYTE );
    memset( x_o2 , 0 , _O2_BYTE );
    memset( temp_o , 0 , sizeof(temp_o) );

    // return: copy w and salt to the signature.
    if( MAX_ATTEMPT_FRMAT <= n_attempt ) return -1;
    memcpy(signature, w, _PUB_N_BYTE);
    memcpy(signature+_PUB_N_BYTE, salt, _SALT_BYTE);
    return 0;
}


static
int _rainbow_verify( const uint8_t * digest , const uint8_t * salt , const unsigned char * digest_ck )
{
    unsigned char correct[_PUB_M_BYTE];
    unsigned char digest_salt[_HASH_LEN + _SALT_BYTE];
    memcpy( digest_salt , digest , _HASH_LEN );
    memcpy( digest_salt+_HASH_LEN , salt , _SALT_BYTE );
    hash_msg( correct , _PUB_M_BYTE , digest_salt , _HASH_LEN+_SALT_BYTE );  // H( digest || salt )

    // check consistancy.
    unsigned char cc = 0;
    for(unsigned i=0;i<_PUB_M_BYTE;i++) {
        cc |= (digest_ck[i]^correct[i]);
    }
    return (0==cc)? 0: -1;
}


int rainbow_verify( const uint8_t * digest , const uint8_t * signature , const pk_t * pk )
{
    unsigned char digest_ck[_PUB_M_BYTE];
    // public_map( digest_ck , pk , signature ); Evaluating the quadratic public polynomials.
    batch_quad_trimat_eval( digest_ck , pk->pk , signature , _PUB_N , _PUB_M_BYTE );

    return _rainbow_verify( digest , signature+_PUB_N_BYTE , digest_ck );
}



///////////////  cyclic version  ///////////////////////////


int rainbow_sign_cyclic( uint8_t * signature , const csk_t * csk , const uint8_t * digest )
{
    sk_t _sk;
    sk_t * sk = &_sk;
    generate_secretkey_cyclic( sk, csk->pk_seed , csk->sk_seed );   // generating classic secret key.

    int r = rainbow_sign( signature , sk , digest );
    memset( sk , 0 , sizeof(sk_t) );  // clean
    return r;
}

int rainbow_verify_cyclic( const uint8_t * digest , const uint8_t * signature , const cpk_t * _pk )
{
    unsigned char digest_ck[_PUB_M_BYTE];
#if defined(_USE_MEMORY_SAVE_)
    // public_map( digest_ck , pk , signature ); Evaluating the quadratic public polynomials.
    rainbow_evaluate_cpk( digest_ck , _pk , signature ); // use less temporary space.
#else
    pk_t _buffer_pk;
    pk_t * pk = &_buffer_pk;
    cpk_to_pk( pk , _pk );         // generating classic public key.
    batch_quad_trimat_eval( digest_ck , pk->pk , signature , _PUB_N , _PUB_M_BYTE );
#endif

    return _rainbow_verify( digest , signature+_PUB_N_BYTE , digest_ck );
}


