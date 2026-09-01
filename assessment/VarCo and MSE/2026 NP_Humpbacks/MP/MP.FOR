      IMPLICIT NONE
      
      INTEGER Yr1, Yr2
      INTEGER Ifeed
      REAL*8 Need
      REAL*8 Catches(1970:2200)
      INTEGER NEstimate
      REAL*8 SurveyY(1000),SurveyEst(1000),SurveyCV(1000)
C
C     Local variables
      INTEGER II,IYR,JJ,KK
      REAL*8 Y1,Y2,StrikeLimit,FinalEst
      
      OPEN(UNIT=13,FILE="MP.Inp")
      READ(13,*)
      READ(13,*) Yr1,Yr2
      READ(13,*)
      READ(13,*) Ifeed
      READ(13,*)
      READ(13,*) Need
      READ(13,*)
      DO 10000 II=Yr1,Yr2
       READ(13,*) IYR,Catches(II)
10000 CONTINUE
      READ(13,*)
      READ(13,*)
      READ(13,*)
      READ(13,*) Nestimate
      READ(13,*)
      DO 20000 II = 1,Nestimate
       READ(13,*) Y1,Y2,(KK,JJ=1,5),SurveyEst(II),SurveyCV(II)
       SurveyY(II) = (Y1+Y2)/2.0
20000 CONTINUE
      CLOSE(13)
      
      FinalEst = SurveyEst(Nestimate)
      IF (FInalEst*0.05.GE.Need) THEN
       StrikeLimit = Need
      ELSE 
       StrikeLimit = FInalEst*0.05
      ENDIF
      WRITE(*,*) Need,FInalEst,StrikeLimit
      
      OPEN(UNIT=14,FILE="MP.OUT")
      WRITE(14,*) StrikeLimit
      WRITE(14,*) Need,FinalEst
      CLOSE(14)

      STOP
      END
      