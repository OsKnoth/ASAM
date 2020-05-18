#!/bin/sh
 set -e # Erwirkt Abbruch des Skripts sobald ein Fehler auftritt

# Aufpassen manche Ordner werden vor der Neuerstellung gelöscht, dies betrifft allerdings nur solche,
# die mit diesem Skript erzeugt werden, daher ist lediglich auf Namengleichheit zu achten
EXAMPLE=Saulesee # Name des zu kompilierenden Beispiels
EX=Saule # Name der Griddatei, alle Ausgabedateien (Weight, Output) sollten den gleichen Namen haben
LAUF=Saulesee_work # Verzeichnis in dem der Modelllauf erfolgte, häufig der gleiche wie das EXAMPLE, in dem Fall:
#LAUF=$EXAMPLE
SRC=SRC_V2.6pre_error # Verzeichnis in dem der ASAM-Quellcode steht
TIME=hou # sek=Sekunde min=Minute hou=Stunde day=Tage # genwünschte Zeiteinheit (nur für Zeitreihe)
READPROF=PROFREAD # Ordner, in dem die Profildateien abgelegt werden
READZEIT=ZEITREAD # Ordner, in dem die Zeitreihen abgelegt werden
UEBERSCHREIBEN="j" # Ordnermanagment (betrifft den Ordner 'ERGEB' :
# 'j': bestehende Ornder überschrieben
# 'n': neue durchnummerierte Ordner erstellen

GRID=$EX.grid
OUT=$EX.out
STARTPFAD=/data/barthel # Auswahl des Ortes, auf dem der ASAM-Quellcode und die Rechnung liegen,
# bei verschiedenen Orten ist dies entsprechen anzupassen
#PFAD=$STARTPFAD/ParExample # Pfad in dem die Rechnung erfolgte meistens in ParExample
PFAD=/home/barthel/Rechnung # Pfad in dem die Rechnung erfolgte meistens in ParExample
ASAM=$STARTPFAD/ASAM # Pfad in dem der ASAM-Quellcode steht
AUS=$PFAD/AUSGABE # Verzeichnis in dem Dateien landen, die die Konsolenausgaben enthalten
MODELLAUF=$PFAD/$LAUF
DATEI=$MODELLAUF/MetOutput # Pfad und Name der Datei, die die Namen der auszugebenden Spezies enthält
ERGEB=Met # Benennt einen Unterordner in dem die Profile nach Namen sortiert abgespeichert werden
ERGEBNIS=$MODELLAUF/$ERGEB  # weist dem entsprechenden Unterordner einem Pfad zu

if [ ! -d $AUS ]; then
  mkdir $AUS
fi

# Kompilieren von ReadProf aus ASAM
cd $ASAM/$SRC
echo Make ReadProf
make ReadProf MACH=gf EX=$EXAMPLE > $AUS/makeProf.out

# nur notwendig, wenn nicht im Exampleordner gerechnet wird
if [ $MODELLAUF != $PFAD/$EXAMPLE ] ; then
  cp -p $PFAD/$EXAMPLE/ReadProf.gf $MODELLAUF
fi

# Erzeugen des Ausgabeverzeichnisses
if [ $UEBERSCHREIBEN = "j" ]; then
  if [ -d $ERGEBNIS ]; then
    rm -r $ERGEBNIS
  fi
else
  j=9999999999
  if [ -d $ERGEBNIS ]; then
    i=0
    while [ $i -lt $j ]
    do
      i=`expr $i + 1`
      if [ ! -d $ERGEBNIS\_$i ]; then
        ERGEBNIS=$ERGEBNIS\_$i
        break
      fi
    done
  fi
fi
mkdir $ERGEBNIS

# In dieser Schleife werden für alle Spezies deren Namen in der oben aufgeführten Datei enthalten sind die Profile und Zeitreihen erzeugt
# Sollte die Schleife nicht verwendet werden, muss sie auskommentiert und der VARIABLE direkt der gewünschte Name zugeordnet werden
# Hierbei werden für chemische Ausgaben die Namen der Spezies wie in der Outputdatei angegeben
# Die meteorologischen Größen wie Temperatur, Feuchte, Wind, usw werden mit der Bezeichnung 'Met' ausgewählt
# Für die SpecialOutputs gilt die Bezeichnung: 'Special'
# Sollten eigene weitergehende Berechnung erfolgen, so ist dies in ReadOutputProfile in der entsprechenden ELSE IF - Schleife möglich
# Hierfür ist die folgende Bezeichnung anzugeben: 'Besonderes'
# Die Nutzung anderer Variablenbezeichnungen führt zu Fehlern

while read LINE
do
  VARIABLE=$LINE
  echo $VARIABLE
  PROF=$ERGEBNIS/$VARIABLE/$READPROF
  ZEIT=$ERGEBNIS/$VARIABLE/$READZEIT
# die nächsten 3 Ordner werden nur im Falle der Flüssigphase benötigt und auch nur dann erzeugt
  MODE=$ERGEBNIS/$VARIABLE/tzC_plot
  CELL=$ERGEBNIS/$VARIABLE/trC_plot
  TIMEORDNER=$ERGEBNIS/$VARIABLE/rzC_plot

#Erzeugen der Verzeichnisse für die Profile und Zeitreihen
  mkdir $ERGEBNIS/$VARIABLE
  mkdir $PROF
  mkdir $ZEIT

  echo 'profile' $VARIABLE

  cd $MODELLAUF

# Zählen der Output-Dateien, die das Modell erzeugt hat
# Notwendig, damit das Profil für alle Dateien erzeugt wird
  Anz=`ls $OUT* | wc -l`
  i=1000
  j=`expr $i + $Anz - 2`

#Erstellen des Profils für jede einzelne Output-Datei
  while [ $i -le $j ]
  do
    cp -p $OUT$i Input.out
    ./ReadProf.gf $GRID $VARIABLE  > $AUS/prof.out # hier muss in dem Ausführungsbefehl an zweiter Stelle der gewünschte Ausgabeparameter stehen, Bsp: Met
    mv Gnu_ProfileMid $PROF/Gnu_Profile$i # schiebt die erzeugten Profile in den vorgesehenen Unterordner
    i=`expr $i + 1`
  done

  cd $PROF

  echo 'Zeitreihe' $VARIABLE
  ls Gnu_Profile* > INPUT_ZEITMITTELUNG # Schreiben der Namen der Profildateien in neue Datei; notwendig für die Erzeugung der Zeitreihen 

  cp -p $ASAM/$SRC/ZEITPROF/Zeitreihe.f90 .
# Kopiert f90.-Datei zum Erzeugen der Zeitreihe
  cp -p $MODELLAUF/$GRID .
# Kopiert Griddatei in aktuelles Verzeichnis, da sie für die Zeitreihenerzeugung relevant ist

# Berechnen der Zeitreihen
  gfortran -O3 -o ausfuehr Zeitreihe.f90
  ./ausfuehr $GRID $TIME $VARIABLE > $AUS/zeit.out # Ausführen der Zeitreihenberechnung mit den Zusatzparamtern Name der Griddatei, Einheit für die Zeitangabe, und Name der ausgegebenen Größe

  mv $PROF/Gnu_timeaverage_cell_* $ZEIT # Schiebt die erzeugten Zeitreihen in den dafür vorgesehenen Ordner

  rm ausfuehr Zeitreihe.f90 $GRID INPUT_ZEITMITTELUNG #löschen nicht mehr benötigter Dateien

  cd $PROF/..

# erzeugt für weitere Ausgaben aus den Zeitreihen Ordner und schiebt die Dateien in selbige
# ermöglicht die Gaskonzentration zeitlich und höhenaufgelöst in 3d darzustellen 
# Gnuplot optimiert
  if [ -f $PROF/Gas_3d_Darstellung ]; then
    mv $PROF/Gas_3d_Darstellung .
  fi

# die folgenden ermöglichen die 3d-Darstellung der Flüssigphase in den verschiedenen Auftragungen entsprechend dem 2 Teil der Dateinamen
# Gnuplot optimiert
  if [ -f $PROF/Gnu_modeaverage_tzC_1 ]; then
    mkdir $MODE
    mv $PROF/Gnu_modeaverage_tzC_* $MODE
  fi

  if [ -f $PROF/Gnu_cellaverage_trC_1 ]; then
    mkdir $CELL
    mv $PROF/Gnu_cellaverage_trC_* $CELL
  fi

  if [ -f $PROF/Gnu_timeaverage_rzC_1 ]; then
    mkdir $TIMEORDNER
    mv $PROF/Gnu_timeaverage_rzC_* $TIMEORDNER
  fi

  cd $MODELLAUF

  rm Input.out # Entfernen nicht mehr benötigter Dateien

  echo 'fertig' $VARIABLE
# nächste Zeile ist im Falle des Nichstverwendens der Schleife auszukommentieren
done < $DATEI

echo 'fertig'
