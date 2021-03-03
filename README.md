# Michael-pakken
Til alle dem der savner Gym-pakken.

# Installation
Læg `MichaelPakken.mla` og evt. også `MichaelPakken.help` ind i Maples lib-mappe.

Stien til lib-mappen kan findes ved at skrive `libname` i en matboks i Maple.
På Windows kan den f.eks. hedde `C:\Program Files\Maple 2020\lib`.
På Mac OS X kan den efter sigende f.eks. hedde `/Library/Frameworks/Maple.framework/Versions/2018/lib`.

`MichaelPakken.mla` gør det muligt at kalde almindeligt brugte funktioner, se en oversigt nedenfor. 
`MichaelPakken.help` gør det muligt at søge igennem alle Mapledemoerne med Maples søgefunktion eller ved at skrive `?søgeord` i Maple. Der ligger også dokumentation for nogle af funktionerne i Michael-pakken.

# Load pakken
Pakken loades ved at skrive `with(Michael)` i en matboks.
Hvis du ønsker at pakken skal loades automatisk hver gang Maple startes (og hver gang der køres en *restart*) kan du tilføje `with(Michael)` til Maples .ini-fil.
På Windows er den fulde sti til .ini-filen `C:\Program Files\Maple 2020\Users\maple.ini`. Den kan oprettes hvis den ikke findes i forvejen. Filen kan redigeres i Notepad el.lign..

# Funktioner i pakken
- **prik(x,y)**: Udregn prikprodukt mellem vektorer (virker for datatyperne liste, vektor og mængde).
- **kryds(x,y)**: Udregn krydsproduktet mellem 3D-vektorer (virker for datatyperne liste, vektor og mængde).
- **længde(x)**: Beregner alm. Euklidisk længde af en vektor (virker for datatyperne liste, vektor og mængde).
- **vop(x)**: Udtag elementerne fra en vektor (virker for datatyperne liste, vektor og mængde).
- **det(A)**: Kort for LinearAlgebra[Determinant](A).
- **gradient(expr)**: Beregner gradienten. Forsøger selv at bestemme variablene i expr, men hvis det ikke virker kan man tilføje en liste, der fortæller hvilke variable, der er tale om, som andet argument.
- **Hessematrix(expr)**: Beregner Hessematricen. Forsøger selv at bestemme variablene i expr, men hvis det ikke virker kan man tilføje en liste, der fortæller hvilke variable, der er tale om, som andet argument.



# Til udviklere (mig)
Når nye funktioner er skrevet gøres følgende
1. Nye funktioner tilføjes som tekst i 'MichaelPakken.mpl' i modulet 'Michael'. Navnet på den nye funktion tilføjes til linjen `export` i toppen.
2. Filen CreateMichaelPakke køres fra start til slut (den flytter den gamle MichaelPakke til mappen Legacy og gemmer en ny 'MichaelPakken.mla'-fil i hoved-mappen).
3. Ændringerne Comittes med Git.

# Kommende funktioner
'RowOperationator' - laver flere RowOperations på én gang med simpel notation.
'Jacobian' - finder selv ud af hvor mange dimensioner. Indbyg mulighed for antagelse om reel, større end værdi, mindre end værdi...
'FindQ' - laver Eigenvectors, ortonormaliserer med Gram-Scmidt, og samler i Q med det(Q)=1 og samme rækkefølge hver gang (sorteret efter egenværdier måske? Evt. lav også FindLambda som så giver samme rækkefølge, og OrtDiag, der giver begge)
'vektorligning' Solve for alle koordinater på 1 gang

Gympakken:
Trigonometriske funktioner med grader
'det' - kort determinant
