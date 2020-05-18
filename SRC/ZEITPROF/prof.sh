#!/bin/sh

# Aufpassen manche Ordner werden vor der Neuerstellung gelöscht, dies betrifft allerdings nur solche,
# die mit diesem Skript erzeugt werden, daher ist lediglich auf Namengleichheit zu achten
EXAMPLE=Saulesee # Name des zu kompilierenden Beispiels
EX=Saule # Name der Griddatei, alle Ausgabedateien (Weight, Output) sollten den gleichen Namen haben
LAUF=Voll # Verzeichnis in dem der Modelllauf erfolgte, häufig der gleiche wie das EXAMPLE, in dem Fall:
# LAUF=$EXAMPLE
SRC=SRC_V2.5 # Verzeichnis in dem der ASAM-Quellcode steht
TIME=hou # sek=Sekunde min=Minute hou=Stunde day=Tage # genwünschte Zeiteinheit (nur für Zeitreihe)
READPROF=PROFREAD # Ordner, in dem die Profildateien abgelegt werden

GRID=$EX.grid
OUT=$EX.out
STARTPFAD=/u1/barthel # Auswahl des Ortes, auf dem der ASAM-Quellcode und die Rechnung liegen,
# bei verschiedenen Orten ist dies entsprechen anzupassen
PFAD=$STARTPFAD/ParExample # Pfad in dem die Rechnung erfolgte meistens in ParExample
ASAM=$STARTPFAD/ASAM # Pfad in dem der ASAM-Quellcode steht
AUS=$PFAD/AUSGABE # Verzeichnis in dem Dateien landen, die die Konsolenausgaben enthalten
MODELLAUF=$PFAD/$LAUF
DATEI=$MODELLAUF/GasOutput # Pfad und Name der Datei, die die Namen der auszugebenden Spezies enthält
ERGEB=Gasphase # Benennt einen Unterordner in dem die Profile nach Namen sortiert abgespeichert werden
ERGEBNIS=$MODELLAUF/$ERGEB  # weist dem entsprechenden Unterordner einem Pfad zu

# Kompilieren von ReadProf aus ASAM
cd $ASAM/$SRC
echo Make ReadProf
make ReadProf MACH=gf EX=$EXAMPLE > $AUS/makeProf.out

# Erzeugen des Ausgabeverzeichnisses
if [ -d $ERGEBNIS ]; then
  rm -r $ERGEBNIS
fi
mkdir $ERGEBNIS

# In dieser Schleife werden für alle Spezies deren Namen in der oben aufgeführten Datei enthalten sind die Profile und Zeitreihen erzeugt
# Sollte die Schleife nicht verwendet werden, muss sie auskommentiert und der VARIABLE direkt der gewünschte Name zugeordnet werden
# Hierbei werden für chemische Ausgaben die Namen der Spezies wie in der Outputdatei angegeben
# Die meteorologischen Größen wie Temperatur, Feuchte, Wind, usw werden mit der Bezeichnung: Met ausgewählt
# Hierbei ist darauf zu achten, dass die gewünschten Größen in der Routine 'ReadOutputProfile' in der Datei 'INIT/Output_Mod.f90' ausgegeben werden 
# SpecialOutputs sind ebenfalls über die meteorologischen Größen anzuwählen, hierfür muss ebenfalls 'ReadOutputProfile' angepasst werden
# Die Nutzung anderer Variablenbezeichnungen führt zu Fehlern
while read LINE; do
   VARIABLE=$LINE
echo $VARIABLE
PROF=$ERGEBNIS/$VARIABLE/$READPROF

#Erzeugen der Verzeichnisse für die Profile und Zeitreihen
mkdir $ERGEBNIS/$VARIABLE
mkdir $PROF

echo 'profile' $VARIABLE

cd $MODELLAUF

# nur notwendig, wenn nicht im Exampleordner gerechnet wird
cp $PFAD/$EXAMPLE/ReadProf.gf .

# Zählen der Output-Dateien, die das Modell erzeugt hat
# Notwendig, damit das Profil für alle Dateien erzeugt wird
Anz=`ls $OUT* | wc -l`
i=1000
j=`expr $i + $Anz - 2`

#Erstellen des Profils für jede einzelne Output-Datei
while [ $i -le $j ]
do
cp -p $OUT$i Input.out
ReadProf.gf $GRID $VARIABLE  > $AUS/prof.out # hier muss in dem Ausführungsbefehl an zweiter Stelle der gewünschte Ausgabeparameter stehen, Bsp: Met
mv Gnu_ProfileMid $PROF/Gnu_Profile$i # schiebt die erzeugten Profile in den vorgesehenen Unterordner
  i=`expr $i + 1`
done

rm Input.out ProfileMid # Entfernen nicht mehr benötigter Dateien

echo 'fertig' $VARIABLE
# nächste Zeile ist im Falle des Nichstverwendens der Schleife auszukommentieren
done < $DATEI

echo 'fertig'
