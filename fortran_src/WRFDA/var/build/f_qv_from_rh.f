












     FUNCTION f_qv_from_rh (RH, T_K, RH0, T_K0, P_PA0) RESULT (QV_KG)



















































      IMPLICIT NONE

      REAL T_K , RH , QV_KG
      REAL P_PA0, T_K0, RH0
      
      REAL ES , QS
      REAL P_MB0, ES0, QS0, QV_KG0




      P_MB0 = P_PA0 / 100.



      ES  = 6.112 * 17.67 * 243.5 * T_K /                       &
                    ((T_K0-273.15+243.5)*(T_K0-273.15+243.5)) * &
                    EXP (17.67*(T_K0-273.15)/(T_K0-273.15+243.5))
      ES0 = 6.112 * EXP (17.67*(T_K0-273.15)/(T_K0-273.15+243.5)) 



      QS  = 0.622 * (P_MB0 * ES ) /  &
                   ((P_MB0 - ES0) * (P_MB0 - ES0))
      QS0 = 0.622 * ES0 /(P_MB0-ES0)            



      QV_KG  = 0.01 * (RH0 * QS + RH * QS0)
      QV_KG0 = 0.01 * RH0 * QS0



      IF (QV_KG0 < 0.) QV_KG = 0.

      RETURN

      END FUNCTION f_qv_from_rh


