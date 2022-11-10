PROGRAM MainProg

  USE Parameter_Mod
  USE IOControl_Mod
  USE Parametric_Mod
  USE GridNeu_Mod
  USE OutputWeightBlk_Mod
! USE OutputWeightNullBlk_Mod
  USE OutputOutGmvGNeu_Mod
  USE OutputUnstructuredGrid_Mod
  USE Emission_Mod
  USE Masking_Mod
  USE Tree_Mod
! USE InitTropoOroStruct_Mod
! USE OutputParaView_Mod
! USE OutUpperBoundsLayer_Mod

  IMPLICIT NONE

  !Local Variable
  CHARACTER(80)  :: InputFileName
  CHARACTER(80)  :: FktInputFileName
  CHARACTER(80)  :: OutputFile
  CHARACTER(100) :: ProgName
  INTEGER :: i,j,l,ib,in
  INTEGER :: dealloc_stat2
  TYPE(Box_T), POINTER :: Box

  version="(Version: V 8.9.4.6)"
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
  !Protokoll-File for AnalysisList
  OPEN(UNIT=OutUnitChL,FILE='.'//TRIM(OutputFile)//'.chl',STATUS='unknown')


  CALL InputControl(InputFileName)
  CALL InitParametric
 
  !Parameter Init
  CALL ComputeParameter

  ! Input 
  CALL InputModelBCVel(InputFileName)
! CALL Display_GridFile(InputFileName)
  CALL ReadGrid(InputFileName)
  CALL Set_Domain_OutView
  CALL GetNeighbor
  CALL NeighborBounds
  CALL Allocate(Floor)

  ! Weights 
! CALL Display_ProcWeights(InputFileName)
  CALL ReadWeights(InputFileName)
  CALL InputFunction(InputFileName)
  CALL InputEmission(InputFileName)
  WRITE(*,*) 'NameFunction  ',NameFunction
  SELECT CASE(NameFunction)
    CASE ('Haus','Haus1') 
      WRITE(*,*) 'Vor WriteHausVTK'
      CALL WriteHausVTK(OutputFile)
      WRITE(*,*) 'Nach WriteHausVTK'
      DO i=1,SIZE(Haus)
        IF (InsideBoxHaus(Root,Haus(i))) THEN
          CALL FindBoxPoint(Root,Box,Haus(i)%PM)
          WRITE(*,*) i,'Haus FindBox',Box%Index,Box%TypeCut,Box%Depth
          WRITE(*,*) 'P0',Box%P0
          WRITE(*,*) 'P0',Haus(i)%P0
          WRITE(*,*) 'PM',Haus(i)%PM
          WRITE(*,*) 'P1',Haus(i)%P1
          WRITE(*,*) 'P1',Box%P1
          S2:DO
            IF (ASSOCIATED(Box)) THEN
              IF (ASSOCIATED(Box%List)) THEN
                S1:DO j=1,SIZE(Box%List)
                  IF (Box%List(j)==i) THEN
                    WRITE(*,*) 'Haus in Box'
                    EXIT S2
                  END IF
                END DO S1
              END IF
              Box=>Box%Parent
            ELSE
              EXIT
            END IF
          END DO S2
        END IF  
      END DO  
  END SELECT    

  ! Vertices
! Call Display_InitAllVertices
  CALL InitAllVertices
! Call Display_AnalyzeAllVertices
  WRITE(*,*) 'AnalyzeAllVertices'
  CALL AnalyzeAllVertices

  ! Edges
! CALL Display_InitAllEdges
  !CALL ProtInitAllEdges
  WRITE(*,*) 'InitAllEdges'
  CALL InitAll
  !CALL ProtInitAllEdges
! CALL Display_AnalyzeAllEdges
  WRITE(*,*) 'AnalyzeAllEdges'
  CALL AnalyzeAllEdges

  ! assort Vertex 
  CALL SortVertex
  WRITE(*,*) 'SortVertex'
  CALL SortVertexView
  WRITE(*,*) 'SortVertexView'
! CALL SortVertexOut2
! CALL SortVertexCutPlane
  WRITE(*,*) 'vor SortVertexSoil'
  CALL SortVertexSoil
  WRITE(*,*) 'SortVertexSoil'
  CALL SortVertexOro

  ! Faces
! CALL Display_AnalyzeAllFaces
  WRITE(*,*) 'AnalyzeAllFaces'
  CALL AnalyzeAllFaces

! ! Cells
! CALL Display_InitAllCells
  WRITE(*,*) 'InitAllCells'
  CALL InitAllCells
! CALL Display_AnalyzeAllCells
  WRITE(*,*) 'AnalyzeAllCells'
  CALL AnalyzeAllCells
! CALL Display_AnalyzeAllCellsRand
! CALL AnalyzeAllCellsRand

! !Test-Check 
! !V_8.9.3.8.1.5, Kinzig_surf_tww97_opt1s200.grid
! !CALL WriteCell(Floor(1)%Cell(32,41,12)%Cell,32,41,12) 
!
! ! Write Output Weight --> Weight2 --->ASAM
! CALL Display_OutWeightBlk(OutputFile)
  WRITE(*,*) 'WriteWeightBlk'
  CALL WriteWeightBlk(OutputFile)
  CALL WriteWeightNetCDF(OutputFile)
! IF (WNull) THEN
!   CALL Display_OutWeightNullBlk(OutputFile)
!   CALL WriteWeightNullBlk(OutputFile)
! END IF
! ! Write Ouput Topography ---> GMV
! CALL Set_ViewOutArea

 WRITE(*,*) 'OutGMV',out_type
 IF (out_type=='a') THEN
!    CALL Display_WriteOutputTopo(out_type)
!    OutputType='ascii'
!    !--------------------------------------------
!    CALL WriteAllCellsAsciiGMV(OutputFile)
!    !--------------------------------------------
!    ! Tropo,Orography,Bounds ---> ParaView,ASAM
!    !--------------------------------------------
!    IF(Ptropo.OR.Pgall.OR.Bound.OR.Pbound) THEN
!       CALL InitTropoBoundOroStruct(OutputFile)
!    END IF
!    IF(Ptropo) CALL WriteTropoCellsParaV(OutputFile)
!    IF(Pgall)  CALL WriteTropoOroCellsParaV(OutputFile)
!    IF(Bound .OR. Pbound) THEN
!       CALL WriteUpperBoundsLayer(OutputFile)
!    END IF
!    !---------------------------------------------
!    IF(GCut)   CALL WriteAllCellsCutPlaneAsciiGMV(OutputFile)
!    IF(GCut2)  CALL WriteAllCellsCutPlaneAsciiGMVTest2(OutputFile)
!    !---------------------------------------------
!    IF(GSoil)  CALL WriteAllCellsSoilAsciiGMV(OutputFile)
!    !---------------------------------------------
!    IF(GONull) CALL WriteAllCellsOroVolNullAsciiGMV(OutputFile)
!    IF(GOro)   Call SortVertAllFacesIn
!    IF(GOro)   Call SortVertCutAllCells
!    IF(GOro)   CALL WriteAllCellsOroAsciiGMV(OutputFile)
 ELSE
     WRITE(*,*) 'Binary '
!    CALL Display_WriteOutputTopo(out_type)
     OutputType='binary'
!    !----------------------------------------------
     CALL WriteAllCellsBinGMV(OutputFile)
     WRITE(*,*) 'Output gmv',gmv
     IF (gmv) THEN
       CALL WriteAllCellsBinary(OutputFile)
     END IF  
     IF (vtk) CALL WriteAllCellsVTK(OutputFile)
     IF (vtk) THEN
       WRITE(*,*) 'vor WriteUnresolvedCellsVTK'
      CALL WriteAllCellsVTK(OutputFile)
      CALL WriteAllBlocksVTK(OutputFile)
      CALL WriteUnresolvedCellsVTK(OutputFile)
       WRITE(*,*) 'nach WriteUnresolvedCellsVTK'
     END IF
!    !--------------------------------------------
!    ! Tropo-,Bound,Orography --->ParaView,ASAM
!    !--------------------------------------------
!    IF(Ptropo.OR.Pgall.OR.Bound.OR.Pbound) THEN
!       CALL InitTropoBoundOroStruct(OutputFile)
!    END IF
!    IF(Ptropo) CALL WriteTropoCellsParaV(OutputFile)
!    IF(Pgall)  CALL WriteTropoOroCellsParaV(OutputFile)
!    IF(Bound .OR. Pbound) THEN
!       CALL WriteUpperBoundsLayer(OutputFile)
!    END IF
!    !----------------------------------------------
     WRITE(*,*) 'GCut  ',GCut,OutputFile
     IF(GCut) THEN
       WRITE(*,*) 'WriteAllCellsCutPlaneBinGMV'
       CALL WriteAllCellsCutPlaneBinGMV(OutputFile)
     END IF  
!    !----------------------------------------------
     IF(GSoil)  CALL WriteAllCellsSoilBinGMV(OutputFile)
!    !---------------------------------------------
!    IF(GONull) CALL WriteAllCellsOroVolNullBinGMV(OutputFile)
     IF (GOro) THEN 
       Call SortVertAllFacesIn
       Call SortVertCutAllCells
       CALL WriteAllCellsOroBinGMV(OutputFile)
     END IF  
  END IF
  CALL WriteUnstructuredGrid(OutputFile)

  IF (Emission) THEN
    CALL ReadEmission
    CALL AnalyzeEmission
    CALL WriteEmission(OutputFile)
  END IF  
  ! Calculate and Write Mask (Baum)
  IF(nr_wahlmask>0) THEN
    Write(*,*) "AnalyzeMasking..."
    CALL AnalyzeMasking !dh
    Write(*,*) "WriteMasking..."
    CALL WriteMasking(OutputFile)
  END IF
 STOP 'nach SortVertAllFacesIn = fertig!' !Roxy, 23.4.2018

! !Domain Checkups
! !Call OutputCheckCellen

! !Write(*,*) "LandDef:"
! !DO i=1,nr_landdef
! !   Write(*,*) "i=",i,"  ",LandDef(1,1,i),LandDef(2,1,i),LandDef(1,2,i),LandDef(2,2,i),LandDef(1,3,i)
! !END DO 

! !Deallocation
! CALL Display_Deallocation
! DO ib=1,nb
!   CALL DeAllocate(Floor(ib))
! END DO  !ib
! Deallocate(Floor)
! !DEALLOCATE VertOut_Lists
! !
! Call DEALLOC_Var_Sierra(nr_wahlfkt)
! Call DEALLOC_Var_SurfitProx(nr_wahlfkt)
! !
! if( nr_landdef>0) then
!   DEALLOCATE(LandDef, STAT=dealloc_stat2)
! end if
! ! Abfrage Mac und ausblenden von Prot_ProgTime()
! ! Ausblenden Prot_ProgTime() -unter Mac-Anwendung
! !                            -unter Linux, wenn gfortran-Debugger
! ! Funktion okey, AIX
! ! CALL Prot_ProgTime()
! CLOSE(UNIT=OutUnitProt)
! CLOSE(UNIT=OutUnitChL)
! CALL Display_InfoProtokollAnalysis(OutputFile)
! CALL Display_InfoEndeProg(version)

END PROGRAM MainProg
