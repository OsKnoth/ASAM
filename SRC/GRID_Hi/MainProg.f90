PROGRAM MainProg

  USE Parameter_Mod
  USE IOControl_Mod
  USE Parametric_Mod
  USE Grid_Mod
  USE OutputWeightBlk_Mod
  USE OutputOutGmvG_Mod

  IMPLICIT NONE

  !Local Variable
  CHARACTER(80)  :: InputFileName
  CHARACTER(80)  :: FktInputFileName
  CHARACTER(80)  :: OutputFile
  CHARACTER(100) :: ProgName
  !CHARACTER(30)  :: version
  INTEGER :: i,l,ib
  INTEGER :: dealloc_stat2              ! return status deallocate

! ========= Input (from READSHAPE) --> Output (for ASAM) ===========
!
!                        *.Haus    -->    *.Weight , ...
!                        *.Street  -->    *.Emission
!                        *.Baum    -->    *.Mask
!
! ACHTUNG! Input-File enthaelt jetzt Objektnummern
!   --> nachtraegliches Eintragen von Minuszeichen vor den Nummern
!       verhindert die Analyse der betreffenden Objekte
!
! ==================================================================

  version="(Version: V 8.9.3.8.1.8)"
  !............................................................................
  !Program Argumente
  !............................................................................
  CALL getarg(0,ProgName)
  CALL getarg(1,InputFileName)
  !.................................
  OutputFile=InputFileName
  i=SCAN(OutputFile,'.',BACK=.TRUE.)
  IF(i>0) THEN
    l=LEN(OutputFile)
    OutputFile(i:l)=" "
  END IF
  !.................................
  !Display Info-Program
  CALL Display_Progname_and_InputFile(ProgName,InputFileName)
  CALL Display_InfoProg(version)

  !Protokoll-File for Analysis, Weight- and GMV-Output 
  OPEN(UNIT=OutUnitProt,FILE=TRIM(OutputFile)//'.prot',STATUS='unknown')
  CALL Prot_Progname_and_InputFile(ProgName,InputFileName)
  CALL Prot_InfoProg(version)

  !Parameter Init
  CALL ComputeParameter

  ! Input 
  CALL Display_GridFile(InputFileName)
  CALL ReadGrid(InputFileName)
  CALL Set_Domain_OutView
  CALL get_ngbrs
  CALL ngbr_bounds
  CALL Allocate(Floor)

  ! Weights 
  CALL Display_ProcWeights(InputFileName) 
  CALL ReadWeights(InputFileName)
  CALL SetWeightBorders

  ! Vertices
  Call Display_InitAllVertices
  CALL InitAllVertices
  Call Display_AnalyzeAllVertices
  CALL AnalyzeAllVertices

  ! Edges
  CALL Display_InitAllEdges
  CALL InitAllEdges
  CALL Display_AnalyzeAllEdges
  CALL AnalyzeAllEdges

  ! assort Vertex 
  CALL SortVertex
  CALL SortVertexCutPlane
  CALL SortVertexSoil
  CALL SortVertexOro

  ! Faces
  !CALL InitAllFaces
  CALL Display_AnalyzeAllFaces
  CALL AnalyzeAllFaces

  ! Cells
  CALL Display_InitAllCells
  CALL InitAllCells
  CALL Display_AnalyzeAllCells
  CALL AnalyzeAllCells
  CALL Display_AnalyzeAllRandCells
  CALL AnalyzeAllRandCells

  !Test-Check 
  !V_8.9.3.8.1.5, Kinzig_surf_tww97_opt1s200.grid
  !CALL WriteCell(Floor(1)%Cell(32,41,12)%Cell,32,41,12) 

  ! Calculate Emission
  IF(nr_wahlemi>0) THEN
    Write(*,*) "AnalyzeEmission..."
    CALL AnalyzeEmission
  END IF

  ! Write Output Weight --> Weight2
  CALL Display_OutWeightBlk(OutputFile)
  CALL WriteWeightBlk(OutputFile)

  ! Write Emission
  IF(nr_wahlemi>0) THEN
    Write(*,*) "WriteEmission..."
    CALL WriteEmission(OutputFile)
  END IF

  ! Calculate and Write Mask (Baum)
  IF(nr_wahlmask>0) THEN
    Write(*,*) "AnalyzeMasking..."
    CALL AnalyzeMasking !dh
    Write(*,*) "WriteMasking..."
    CALL WriteMasking(OutputFile)
  END IF

  ! Write Ouput Topography ---> GMV
  CALL Set_ViewOutArea
  IF (out_type=='a') THEN
     CALL Display_WriteOutputTopo(out_type)
     OutputType='ascii'
     !--------------------------------------------
     CALL WriteAllCellsAsciiGMV(OutputFile)
     !---------------------------------------------
     CALL WriteAllCellsCutPlaneAsciiGMV(OutputFile)
     CALL WriteAllCellsCutPlaneAsciiGMVTest2(OutputFile)
     !---------------------------------------------
     CALL WriteAllCellsSoilAsciiGMV(OutputFile)
     !---------------------------------------------
  !  !Call SortVertAllFacesIn
  !  !Call SortVertCutAllCells
  !  !CALL WriteAllCellsOroAsciiGMV(OutputFile)
  ELSE
     CALL Display_WriteOutputTopo(out_type)
     OutputType='binary'
     !----------------------------------------------
     CALL WriteAllCellsBinGMV(OutputFile)
     !----------------------------------------------
     CALL WriteAllCellsCutPlaneBinGMV(OutputFile)
     !----------------------------------------------
     CALL WriteAllCellsSoilBinGMV(OutputFile)
     !---------------------------------------------
  !   Call SortVertAllFacesIn
  !   Call SortVertCutAllCells
  !   CALL WriteAllCellsOroBinGMV(OutputFile)
  END IF

  !Domain Checkups
  !Call OutputCheckCellen

  !Write(*,*) "LandDef:"
  !DO i=1,nr_landdef
  !   Write(*,*) "i=",i,"  ",LandDef(1,1,i),LandDef(2,1,i),LandDef(1,2,i),LandDef(2,2,i),LandDef(1,3,i)
  !END DO 

  !Deallocation
  CALL Display_Deallocation
  DO ib=1,nb
    CALL DeAllocate(Floor(ib))
  END DO  !ib
  Deallocate(Floor)
  !
  Call DEALLOC_Var_Sierra(nr_wahlfkt)
  Call DEALLOC_Var_SurfitProx(nr_wahlfkt)
  !
  if( nr_landdef>0) then
    DEALLOCATE(LandDef, STAT=dealloc_stat2)
  end if

  CALL Prot_ProgTime()
  CLOSE(UNIT=OutUnitProt)
  CALL Display_InfoProtokollAnalysis(OutputFile)
  CALL Display_InfoEndeProg(version)

END PROGRAM MainProg
