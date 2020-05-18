MODULE Sp_Mod

  USE Kind_Mod
  USE Tool_Mod
  USE MatSpDiag_Mod
  USE MatSpRowCol_Mod
  USE MatSpRowColMat_Mod
  USE MatSpRowColD_Mod
  USE MatSpRowColDiag_Mod
  USE MatSpRowColL_Mod
  USE MatSpInd_Mod
  USE MatSpConvert_Mod

  INTERFACE SpNullify
    MODULE PROCEDURE SpNullify_SpRowCol  &
                    ,SpNullify_SpRowColMat  &
                    ,SpNullify_SpDiag    &
                    ,SpNullify_SpRowColL &
                    ,SpNullify_SpRowColDiag &
                    ,SpNullify_SpInd     &
                    ,SpNullify_SpIndMat     &
                    ,SpNullify_SpRowColD 
  END INTERFACE
  INTERFACE SpInsert
    MODULE PROCEDURE Insert_SpRowColD 
  END INTERFACE
  INTERFACE gcmat
    MODULE PROCEDURE gcmat_SpRowColD
  END INTERFACE
  INTERFACE SpAVec
    MODULE PROCEDURE SpAVec_SpRowCol &
                    ,SpAVecI4_SpRowColMat &
                    ,SpAVec_SpDiag   &
                    ,SpABVec_SpDiag   &
                    ,SpAVec1_SpDiag   &
                    ,SpABVec1_SpDiag   &
                    ,SpAVecI3_SpRowColL &
                    ,SpAVecI3_SpRowCol &
                    ,SpAVecI4_SpRowCol &
                    ,SpAVec_SpRowColL
  END INTERFACE
  INTERFACE SymbLU
    MODULE PROCEDURE SymbLU_SpRowColD
  END INTERFACE
  INTERFACE Axpy
    MODULE PROCEDURE Axpy_SpRowColL &
                    ,Axpy_SpRowCol &
                    ,AddDiag3_SpRowCol &
                    ,AddDiag_SpRowCol &
                    ,AddDiag3_SpRowColMat
  END INTERFACE
  INTERFACE ATxpy
    MODULE PROCEDURE ATxpy_SpRowCol
  END INTERFACE
  INTERFACE GaussSeidel
    MODULE PROCEDURE GaussSeidel_SpRowColL &
                    ,GaussSeidel_SpDiag &
                    ,GaussSeidel_SpDiagAB
  END INTERFACE
  INTERFACE GaussSeidelF
    MODULE PROCEDURE GaussSeidelF_SpRowColL
  END INTERFACE
  INTERFACE GaussSeidelB
    MODULE PROCEDURE GaussSeidelB_SpRowColL
  END INTERFACE
  INTERFACE IC0
    MODULE PROCEDURE IC_SpRowColL
  END INTERFACE
  INTERFACE IC0Solve
    MODULE PROCEDURE ICSolve_SpRowColL &
                    ,ICSolveM_SpRowColL &
                    ,ICSolveMI3_SpRowColL &
                    ,GaussSeidel_SpRowCol &
                    ,GaussSeidelI3_SpRowCol &
                    ,GaussSeidelI4_SpRowColMat 
  END INTERFACE
! INTERFACE ASSIGNMENT(=)SpConvert
  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE RowCol_To_Diag &
                    ,Ind_To_RowCol &
                    ,Ind_To_RowColMat &
                    ,RowCol_To_RowColL & 
                    ,RowColD_To_RowColDiag &  
                    ,Copy_MatSpDiag
  END INTERFACE
  INTERFACE SpMm
    MODULE PROCEDURE SpMm_SpRowCol &
                    ,SpMm_SpRowColMat
  END INTERFACE
  INTERFACE SpTrans
    MODULE PROCEDURE  SpTrans_SpRowCol
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE SpAllocate_SpInd &
                    ,SpAllocate_SpIndMat &
                    ,SpAllocate_SpRowColD
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE SpDeallocate_SpInd     &
                    ,SpDeallocate_SpIndMat  &
                    ,SpDeallocate_SpRowCol  &
                    ,SpDeallocate_SpRowColMat  &
                    ,SpDeallocate_SpRowColDiag  &
                    ,SpDeallocate_SpRowColD
  END INTERFACE
  INTERFACE Output
    MODULE PROCEDURE  SpOutput_SpDiag,SpOutput_SpRowColL, &
                      SpOutput_SpRowCol
  END INTERFACE

END MODULE Sp_Mod
  

