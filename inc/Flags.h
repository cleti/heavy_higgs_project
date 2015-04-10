
#ifndef FLAGS_H
#define FLAGS_H

#include <boost/utility.hpp>

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

#define I_FLAGS_ALL     BOOST_BINARY(1 111 111 1)
#define I_FLAGS_B_QCD   BOOST_BINARY(0 000 000 1)
#define I_FLAGS_B_PHI   BOOST_BINARY(0 000 001 0)
#define I_FLAGS_B       BOOST_BINARY(0 000 001 1)
#define I_FLAGS_V       BOOST_BINARY(0 000 010 0)
#define I_FLAGS_D       BOOST_BINARY(0 000 100 0)
#define I_FLAGS_R_GG    BOOST_BINARY(0 001 000 0)
#define I_FLAGS_R_QQ    BOOST_BINARY(0 010 000 0)
#define I_FLAGS_R_QG    BOOST_BINARY(0 100 000 0)
#define I_FLAGS_PLOT    BOOST_BINARY(1 000 000 0)
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

#define SET_FLAG(FLAG,F)  F |=  FLAG
#define USET_FLAG(FLAG,F) F &= ~FLAG
  
// B3   | B2   | B1   | CTV2 | V2   | CTV1 | V1   | 4G   | CTSE |  SE  | Born
#define F_EVAL_ALL           BOOST_BINARY_LU(111 111 111 111 111 111 111 1)

#define F_EVAL_B_ALL         BOOST_BINARY_LU(000 000 000 000 000 000 111 0)
#define F_EVAL_B_PHIxPHI     BOOST_BINARY_LU(000 000 000 000 000 000 001 0)
#define F_EVAL_B_PHIxQCD     BOOST_BINARY_LU(000 000 000 000 000 000 010 0)
#define F_EVAL_B_QCDxQCD     BOOST_BINARY_LU(000 000 000 000 000 000 100 0)	     		        

#define F_EVAL_V_ALL         BOOST_BINARY_LU(000 000 111 111 111 111 000 0)
#define F_EVAL_V_PHI0xQCD1   BOOST_BINARY_LU(000 000 001 111 111 111 000 0)
#define F_EVAL_V_SE          BOOST_BINARY_LU(000 000 000 000 000 001 000 0)
#define F_EVAL_V_CTSE        BOOST_BINARY_LU(000 000 000 000 000 010 000 0)
#define F_EVAL_V_4G          BOOST_BINARY_LU(000 000 000 000 000 100 000 0)
#define F_EVAL_V_V1          BOOST_BINARY_LU(000 000 000 000 001 000 000 0)
#define F_EVAL_V_CTV1        BOOST_BINARY_LU(000 000 000 000 010 000 000 0)
#define F_EVAL_V_V2          BOOST_BINARY_LU(000 000 000 000 100 000 000 0)
#define F_EVAL_V_CTV2        BOOST_BINARY_LU(000 000 000 001 000 000 000 0)
#define F_EVAL_V_B1          BOOST_BINARY_LU(000 000 000 010 000 000 000 0)
#define F_EVAL_V_B2          BOOST_BINARY_LU(000 000 000 100 000 000 000 0)
#define F_EVAL_V_B3          BOOST_BINARY_LU(000 000 001 000 000 000 000 0)
#define F_EVAL_V_PHI1xQCD0   BOOST_BINARY_LU(000 000 010 000 000 000 000 0)
#define F_EVAL_V_PHIxPHI     BOOST_BINARY_LU(000 000 100 000 000 000 000 0)
#define F_EVAL_V_NF          BOOST_BINARY_LU(000 001 000 000 000 000 000 0)
					     
#define F_EVAL_D_ALL         BOOST_BINARY_LU(111 000 000 000 000 000 000 0)
#define F_EVAL_D_DELTA       BOOST_BINARY_LU(001 000 000 000 000 000 000 0)
#define F_EVAL_D_END         BOOST_BINARY_LU(010 000 000 000 000 000 000 0)
#define F_EVAL_D_CONT        BOOST_BINARY_LU(100 000 000 000 000 000 000 0)





#define SET_EVAL_ALL(F)     F |= F_EVAL_ALL  

#define SET_EVAL_B_ALL(F)       F |= F_EVAL_B_ALL  
#define SET_EVAL_B_PHIxPHI(F)   F |= F_EVAL_B_PHIxPHI  
#define SET_EVAL_B_PHIxQCD(F)   F |= F_EVAL_B_PHIxQCD  
#define SET_EVAL_B_QCDxQCD(F)   F |= F_EVAL_B_QCDxQCD  

#define SET_EVAL_V_ALL(F)   F |= F_EVAL_V_ALL  
#define SET_EVAL_V_SE(F)    F |= F_EVAL_V_SE   
#define SET_EVAL_V_CTSE(F)  F |= F_EVAL_V_CTSE 
#define SET_EVAL_V_4G(F)    F |= F_EVAL_V_4G   
#define SET_EVAL_V_V1(F)    F |= F_EVAL_V_V1   
#define SET_EVAL_V_CTV1(F)  F |= F_EVAL_V_CTV1 
#define SET_EVAL_V_V2(F)    F |= F_EVAL_V_V2   
#define SET_EVAL_V_CTV2(F)  F |= F_EVAL_V_CTV2 
#define SET_EVAL_V_B1(F)    F |= F_EVAL_V_B1   
#define SET_EVAL_V_B2(F)    F |= F_EVAL_V_B2   
#define SET_EVAL_V_B3(F)    F |= F_EVAL_V_B3
#define SET_EVAL_V_PHI0xQCD1(F) F |= F_EVAL_V_PHI0xQCD1
#define SET_EVAL_V_PHI1xQCD0(F) F |= F_EVAL_V_PHI1xQCD0
#define SET_EVAL_V_PHIxPHI(F)   F |= F_EVAL_V_PHIxPHI
#define SET_EVAL_V_NF(F)        F |= F_EVAL_V_NF
  
#define SET_EVAL_D_ALL(F)    F |= F_EVAL_D_ALL
#define SET_EVAL_D_DELTA(F)  F |= F_EVAL_D_DELTA
#define SET_EVAL_D_END(F)    F |= F_EVAL_D_END   
#define SET_EVAL_D_CONT(F)   F |= F_EVAL_D_CONT   



#define USET_EVAL_ALL(F)     F &= ~F_EVAL_ALL  

#define USET_EVAL_B_ALL(F)       F &= ~F_EVAL_B_ALL  
#define USET_EVAL_B_PHIxPHI(F)   F &= ~F_EVAL_B_PHIxPHI  
#define USET_EVAL_B_PHIxQCD(F)   F &= ~F_EVAL_B_PHIxQCD  
#define USET_EVAL_B_QCDxQCD(F)   F &= ~F_EVAL_B_QCDxQCD  

#define USET_EVAL_V_ALL(F)   F &= ~F_EVAL_V_ALL  
#define USET_EVAL_V_SE(F)    F &= ~F_EVAL_V_SE   
#define USET_EVAL_V_CTSE(F)  F &= ~F_EVAL_V_CTSE 
#define USET_EVAL_V_4G(F)    F &= ~F_EVAL_V_4G   
#define USET_EVAL_V_V1(F)    F &= ~F_EVAL_V_V1   
#define USET_EVAL_V_CTV1(F)  F &= ~F_EVAL_V_CTV1 
#define USET_EVAL_V_V2(F)    F &= ~F_EVAL_V_V2   
#define USET_EVAL_V_CTV2(F)  F &= ~F_EVAL_V_CTV2 
#define USET_EVAL_V_B1(F)    F &= ~F_EVAL_V_B1   
#define USET_EVAL_V_B2(F)    F &= ~F_EVAL_V_B2   
#define USET_EVAL_V_B3(F)    F &= ~F_EVAL_V_B3
#define USET_EVAL_V_PHI0xQCD1(F) F &= ~F_EVAL_V_PHI0xQCD1
#define USET_EVAL_V_PHI1xQCD0(F) F &= ~F_EVAL_V_PHI1xQCD0
#define USET_EVAL_V_PHIxPHI(F)   F &= ~F_EVAL_V_PHIxPHI
#define USET_EVAL_V_NF(F)        F &= ~F_EVAL_V_NF

#define USET_EVAL_D_ALL(F)    F &= ~F_EVAL_D_ALL
#define USET_EVAL_D_DELTA(F)  F &= ~F_EVAL_D_DELTA
#define USET_EVAL_D_END(F)    F &= ~F_EVAL_D_END 
#define USET_EVAL_D_CONT(F)   F &= ~F_EVAL_D_CONT   

// test if any born or virtual flags are set
#define EVAL_B(F)     (F & (F_EVAL_B_ALL))
#define EVAL_B_PHI(F) (F & (F_EVAL_B_PHIxPHI | F_EVAL_B_PHIxQCD))
#define EVAL_V(F)     (F & (F_EVAL_V_ALL))
// test if any int. dipole flags are set
#define EVAL_ID(F)    (F & F_EVAL_D_ALL)

#define IZ_EVAL_SI_2P(F)  (F & F_EVAL_V_SE) || IZ_EVAL_SI_3P(F)
#define IZ_EVAL_SI_3P(F)  (F & (F_EVAL_V_4G | F_EVAL_V_V1 | F_EVAL_V_V2 | F_EVAL_V_PHI1xQCD0 | F_EVAL_V_PHIxPHI)) || IZ_EVAL_SI_4P(F)
#define IZ_EVAL_SI_4P(F)   F & (F_EVAL_V_B1 | F_EVAL_V_B2 | F_EVAL_V_B3)

////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////000 000 000 001 000 0/
#define F_EVAL_R_ALL          BOOST_BINARY_LU(111 111 111 111 111 1)
#define F_EVAL_R_GG           BOOST_BINARY_LU(000 110 011 111 111 1)
#define F_EVAL_R_QQ           BOOST_BINARY_LU(001 000 100 000 000 0)
#define F_EVAL_R_QG           BOOST_BINARY_LU(010 001 000 000 000 0)
//////////////// QCD_PHI ///////////////////////////////////////////
#define F_EVAL_R_ISR_ISR      BOOST_BINARY_LU(000 000 000 000 000 1)
#define F_EVAL_R_ISR_FSR      BOOST_BINARY_LU(000 000 000 000 001 0)
#define F_EVAL_R_ISR_INT      BOOST_BINARY_LU(000 000 000 000 010 0)
#define F_EVAL_R_FSR_ISR      BOOST_BINARY_LU(000 000 000 000 100 0)
#define F_EVAL_R_FSR_FSR      BOOST_BINARY_LU(000 000 000 001 000 0)
#define F_EVAL_R_FSR_INT      BOOST_BINARY_LU(000 000 000 010 000 0)
#define F_EVAL_R_INT_ISR      BOOST_BINARY_LU(000 000 000 100 000 0)
#define F_EVAL_R_INT_FSR      BOOST_BINARY_LU(000 000 001 000 000 0)
#define F_EVAL_R_INT_INT      BOOST_BINARY_LU(000 000 010 000 000 0)
#define F_EVAL_R_PHIxQCD_QQ   BOOST_BINARY_LU(000 000 100 000 000 0)
#define F_EVAL_R_PHIxQCD_QG   BOOST_BINARY_LU(000 001 000 000 000 0)
#define F_EVAL_R_PHIxPHI_ISR  BOOST_BINARY_LU(000 010 000 000 000 0)
#define F_EVAL_R_PHIxPHI_FSR  BOOST_BINARY_LU(000 100 000 000 000 0)
#define F_EVAL_R_PHIxPHI_QQ   BOOST_BINARY_LU(001 000 000 000 000 0)
#define F_EVAL_R_PHIxPHI_QG   BOOST_BINARY_LU(010 000 000 000 000 0)

// unintegrated dipole labels: g1 g2 t tb
// initial-initial dipoles needed to make ISR-ISR real corrections finite
#define F_EVAL_UID_ES00   F_EVAL_R_ISR_ISR
#define F_EVAL_UID_SE00   F_EVAL_R_ISR_ISR

#define F_EVAL_UID_E0S0   F_EVAL_R_ISR_FSR | F_EVAL_R_FSR_ISR
#define F_EVAL_UID_S0E0   F_EVAL_R_ISR_FSR | F_EVAL_R_FSR_ISR
#define F_EVAL_UID_0ES0   F_EVAL_R_ISR_FSR | F_EVAL_R_FSR_ISR
#define F_EVAL_UID_0E0S   F_EVAL_R_ISR_FSR | F_EVAL_R_FSR_ISR
#define F_EVAL_UID_0S0E   F_EVAL_R_ISR_FSR | F_EVAL_R_FSR_ISR
#define F_EVAL_UID_S00E   F_EVAL_R_ISR_FSR | F_EVAL_R_FSR_ISR
#define F_EVAL_UID_E00S   F_EVAL_R_ISR_FSR | F_EVAL_R_FSR_ISR

#define F_EVAL_UID_00ES   F_EVAL_R_FSR_FSR
#define F_EVAL_UID_00SE   F_EVAL_R_FSR_FSR


#define SET_EVAL_R_ALL(F)      F |= F_EVAL_R_ALL  
#define SET_EVAL_R_ISR_ISR(F)  F |= F_EVAL_R_ISR_ISR
#define SET_EVAL_R_ISR_FSR(F)  F |= F_EVAL_R_ISR_FSR
#define SET_EVAL_R_ISR_INT(F)  F |= F_EVAL_R_ISR_INT
#define SET_EVAL_R_FSR_ISR(F)  F |= F_EVAL_R_FSR_ISR
#define SET_EVAL_R_FSR_FSR(F)  F |= F_EVAL_R_FSR_FSR
#define SET_EVAL_R_FSR_INT(F)  F |= F_EVAL_R_FSR_INT
#define SET_EVAL_R_INT_ISR(F)  F |= F_EVAL_R_INT_ISR
#define SET_EVAL_R_INT_FSR(F)  F |= F_EVAL_R_INT_FSR
#define SET_EVAL_R_INT_INT(F)  F |= F_EVAL_R_INT_INT
#define SET_EVAL_R_PHIxQCD_QQ(F)   F |= F_EVAL_R_PHIxQCD_QQ
#define SET_EVAL_R_PHIxQCD_QG(F)   F |= F_EVAL_R_PHIxQCD_QG 
#define SET_EVAL_R_PHIxPHI_ISR(F)  F |= F_EVAL_R_PHIxPHI_ISR
#define SET_EVAL_R_PHIxPHI_FSR(F)  F |= F_EVAL_R_PHIxPHI_FSR
#define SET_EVAL_R_PHIxPHI_QQ(F)  F |= F_EVAL_R_PHIxPHI_QQ
#define SET_EVAL_R_PHIxPHI_QG(F)  F |= F_EVAL_R_PHIxPHI_QG

#define USET_EVAL_R_ALL(F)      F &= ~F_EVAL_R_ALL  
#define USET_EVAL_R_ISR_ISR(F)  F &= ~F_EVAL_R_ISR_ISR
#define USET_EVAL_R_ISR_FSR(F)  F &= ~F_EVAL_R_ISR_FSR
#define USET_EVAL_R_ISR_INT(F)  F &= ~F_EVAL_R_ISR_INT
#define USET_EVAL_R_FSR_ISR(F)  F &= ~F_EVAL_R_FSR_ISR
#define USET_EVAL_R_FSR_FSR(F)  F &= ~F_EVAL_R_FSR_FSR
#define USET_EVAL_R_FSR_INT(F)  F &= ~F_EVAL_R_FSR_INT
#define USET_EVAL_R_INT_ISR(F)  F &= ~F_EVAL_R_INT_ISR
#define USET_EVAL_R_INT_FSR(F)  F &= ~F_EVAL_R_INT_FSR
#define USET_EVAL_R_INT_INT(F)  F &= ~F_EVAL_R_INT_INT
#define USET_EVAL_R_PHIxQCD_QQ(F)   F &= ~F_EVAL_R_PHIxQCD_QQ
#define USET_EVAL_R_PHIxQCD_QG(F)   F &= ~F_EVAL_R_PHIxQCD_QG 
#define USET_EVAL_R_PHIxPHI_ISR(F)  F &= ~F_EVAL_R_PHIxPHI_ISR
#define USET_EVAL_R_PHIxPHI_FSR(F)  F &= ~F_EVAL_R_PHIxPHI_FSR
#define USET_EVAL_R_PHIxPHI_QQ(F)  F &= ~F_EVAL_R_PHIxPHI_QQ
#define USET_EVAL_R_PHIxPHI_QG(F)  F &= ~F_EVAL_R_PHIxPHI_QG

#define EVAL_R_QG(F) (F & (F_EVAL_R_PHIxPHI_QG | F_EVAL_R_PHIxQCD_QG))
#define EVAL_R_QQ(F) (F & (F_EVAL_R_PHIxPHI_QQ | F_EVAL_R_PHIxQCD_QQ))
#define EVAL_R_QQ_QG(F) (EVAL_R_QG(F) | EVAL_R_QQ(F))
#define EVAL_UID_GG_II(F)  (F & (F_EVAL_R_ISR_ISR | F_EVAL_R_PHIxPHI_ISR))
#define EVAL_UID_GG_FF(F)  (F & (F_EVAL_R_FSR_FSR | F_EVAL_R_PHIxPHI_FSR))
#define EVAL_UID_QG_II(F)  EVAL_R_QG(F)

  
#define IZ_EVAL_UID_INT(F)   F & (F_EVAL_R_ISR_FSR | F_EVAL_R_FSR_ISR )
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////




#endif
