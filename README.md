# Michael-pakken
Til alle dem der savner Gym-pakken.

# Installation
Læg `MichaelPakken.mla` og evt. også `MichaelPakken.help` og `Integrator8.mw`  ind i Maples lib-mappe.

Stien til lib-mappen kan findes ved at skrive `libname` i en matboks i Maple.
På Windows kan den f.eks. hedde `C:\Program Files\Maple 2020\lib`.
På MacOS kan den hedde `/Library/Frameworks/Maple.framework/Versions/2020/lib`.

`MichaelPakken.mla` gør det muligt at kalde almindeligt brugte funktioner, se en oversigt nedenfor.
`MichaelPakken.help` gør det muligt at søge igennem alle Mapledemoerne med Maples søgefunktion eller ved at skrive `?søgeord` i Maple. Der ligger også dokumentation for nogle af funktionerne i Michael-pakken. Hvis man ikke søger på det præcise navn for Mapledemoen er det meget sandsynligt, at nogle andre hjælpesider åbnes op først. I så fald kan man rulle ned til 'User Help' ude i siden af Maples søgevindue, og se hvilke hjælpefiler, der ligger dér.

# Load pakken
Pakken loades ved at skrive `with(Michael)` i en matboks.
Hvis du ønsker at pakken skal loades automatisk hver gang Maple startes (og hver gang der køres en *restart*) kan du tilføje `with(Michael)` til Maples initialiserings-fil.
På Windows er den fulde sti til .ini-filen `C:\Program Files\Maple 2020\Users\maple.ini`.
På MacOS er den `/Library/Frameworks/Maple.framework/Versions/2020/lib/init` (bemærk ingen file extension).
Filen kan oprettes hvis den ikke findes i forvejen. Filen kan redigeres i Notepad el.lign..

# Funktioner i pakken
Næsten alle vektorfunktioner fungerer også, hvis der passes en liste eller mængde i stedet.
- **vop(x)**: Udtager elementerne fra en vektor (en generalisering af **op**).
- **prik(x,y)**: Udregn prikprodukt mellem vektorer.
- **prikc(x,y)**: Udregn det standard hermitiske indre produkt mellem vektorer (dvs. **prik** inkl. konjugering).
- **kryds(x,y)**: Udregn krydsproduktet mellem 3D-vektorer.
- **rum(x,y,z)**: Rumprodukt mellem 3 3D-vektorer.
- **længde(x)**: Beregner alm. Euklidisk længde af en vektor.
- **grad(expr,variabel_liste)**: Beregner gradienten af en vektor med hensyn til angivne variable.
- **det(A)**: Kort for LinearAlgebra\[Determinant\](A).
- **div(V)**: Divergens for vektorfelt.
- **rot(V)**: Rotation for vektorfelt.
- **Hessematrix(expr, variabel_liste)**: Beregner Hessematricen.
- **GetJacobi(parameterfremstilling, variabel_liste)**: Beregner Jacobi for en parameterfremstilling. Finder selv ud af, hvilken slags Jacobi-fkt., der er tale om. Kan med fordel kombineres med `assuming`, f.eks. `simplify(GetJacobi(r(u,v),[u,v])) assuming u>=0, u<=2*Pi, v>=0`.
- **TrappeMetode(V)**: Beregner kurveintegralet af V fra origo til et vilkårligt punkt (x,y,z). Det kan nemt tjekkes om V er et gradientvektorfelt ved at beregne **grad(TrappeMetode(V),variabel_liste)**, og se, om det er det samme som V.



# Til udviklere (mig)
Når nye funktioner er skrevet gøres følgende
1. Nye funktioner tilføjes som tekst i 'MichaelPakken.mpl' i modulet 'Michael'. Navnet på den nye funktion tilføjes til linjen `export` i toppen.
2. Filen CreateMichaelPakke køres fra start til slut (den flytter den gamle MichaelPakke til mappen Legacy og gemmer en ny 'MichaelPakken.mla'-fil i hoved-mappen).
3. Ændringerne comittes med Git.

Nye hjælpefiler laves på følgende måde
1. Skriv hjælpefilen, enten som .mw eller .txt
2. Gem filen i MichaelPakken/Kodning/Help Database
3. Kør anden sektion i CreateHelpDatabase, som tilføjer nye filer til databasefilen MichaelPakken.help
4. Ændringer comittes med Git

For ikke at behøve at kopiere MichaelPakken.mla og MichaelPakken.help til lib-mappen hver gang de opdateres, har jeg tilføjet flg. til min .ini-fil:
`libname := "sti til dette repository", libname:`


# Kommende funktioner
'RowOperationator' - laver flere RowOperations på én gang med simpel notation.
'Jacobian' - finder selv ud af hvor mange dimensioner. Indbyg mulighed for antagelse om reel, større end værdi, mindre end værdi...
'FindQ' - laver Eigenvectors, ortonormaliserer med Gram-Scmidt, og samler i Q med det(Q)=1 og samme rækkefølge hver gang (sorteret efter egenværdier måske? Evt. lav også FindLambda som så giver samme rækkefølge, og OrtDiag, der giver begge)
'vektorligning' Solve for alle koordinater på 1 gang

Gympakken:
Trigonometriske funktioner med grader
'det' - kort determinant
