# Michael-pakken
Til alle dem der savner Gym-pakken.

# Installation
Læg `MichaelPakken.mla` ind i Maples lib-mappe.

Stien til lib-mappen kan findes ved at skrive `libname` i en matboks i Maple.
På Windows kan den f.eks. hedde `C:\Program Files\Maple 2020\lib`.
På Mac OS X kan den efter sigende f.eks. hedde `/Library/Frameworks/Maple.framework/Versions/2018/lib`.

# Load pakken
Pakken loades ved at skrive `with(Michael)` i en matboks.
Hvis du ønsker at pakken skal loades automatisk hver gang Maple startes (og hver gang der køres en *restart*) kan du tilføje `with(Michael)` til Maples .ini-fil.
På Windows er den fulde sti til .ini-filen `C:\Program Files\Maple 2020\Users\maple.ini`. Den kan oprettes hvis den ikke findes i forvejen. Filen kan redigeres i Notepad el.lign..

# Development
Nye funktioner tilføjes i 'MichaelPakken.mpl'. Den gamle 'MichaelPakken.mla' lægges ind i mappen 'Legacy' under et unikt navn. CreateMichaelPakke køres fra start, så den gemmer en ny 'MichaelPakken.mla'-fil i hoved-mappen. Ændringerne Comittes.

# Kommende funktioner
'RowOperationator' - laver flere RowOperations på én gang med simpel notation.
'Jacobian' - finder selv ud af hvor mange dimensioner. Indbyg mulighed for antagelse om reel, større end værdi, mindre end værdi...
'FindQ' - laver Eigenvectors, ortonormaliserer med Gram-Scmidt, og samler i Q med det(Q)=1 og samme rækkefølge hver gang (sorteret efter egenværdier måske? Evt. lav også FindLambda som så giver samme rækkefølge, og OrtDiag, der giver begge)
'vektorligning' Solve for alle koordinater på 1 gang

Gympakken:
Trigonometriske funktioner med grader
'det' - kort determinant
