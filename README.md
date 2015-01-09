projekt-mapper
==============
Sustav za mapiranje DNA sekvenci na referentni genom


<h3>Opis:</h3>

<p>
Ovaj rad opisuje sustav za mapiranje kratkih očitanja na poznati genom.
Potrebno je svako kratko očitanje (read) koje su nalazi u FASTQ datoteci mapirati na referentni genom, 
Kako je prilikom određivanja kratkih očitanja moglo doći do mutacija ili različitih pogrešaka, zadaća je odrediti najizgledniju poziciju zadanog očitanja unutar referentnog genoma. Zbog veličine genoma nije isplativo provođenje lokalnog poravnavanja na čitavom genomu, stoga je cilj odrediti regije u referentnom genomu koje imaju veću vjerojatnost da upravu njima pripada očitanje.
</p>
<p>
Ovaj program to radi na sljedeći način. 
Read se podjeli na sve podnizove duljine k koji se naziva kmer. Kmer je jednoznacno određen pozicijom na kojom započinje unutar read-a.
Za svaki kmer odredimo pozicije pojavljivanja unutar referentnog genoma. Time u konačnici imamo polje parova (indeks početka u referentnom genomu, indeks početka u read-u odnosno identifikator kmera). Polje parova sortira se po prvoj vrijednosti (indeksu u referentnom genomu).<br>
Lako je zaključiti da regije u kojima je vjerojatnije da se to očitanje sadržavati uzastopne kmerove.
Tu se javlja problem odabira adekvatne metrike koja će odrediti kandidatne regije.
Jedna od metrika duljina rastućeg podniza. Pod rastući podniz misli se na rastuće indekse kmer-a. 
Druga metrika je LCSk. LCSk predstavlja broj nepreklapajućih podnizova duljine k između daju nizova znakova. Ova metrika se pokazala u praksi bolja
prilikom određivanja kandidatnih regija.
Kako bi se ograničilo promatrano područje promatraju se oni parovi koji upadaju u prozor dva puta dulji od duljine read-a.<br>
Nakon sta odaberemo regiju s najvećim LCSk, tada se na njoj provede Smith-Waterman algoritam za lokalno poravnavanje.


<h3>Potrebne dodatne biblioteke:</h3>
libdivsufsort - https://code.google.com/p/libdivsufsort/

<hr>Korištenje:</h3>
<code>mapper construct 'fastaFile' 'suffixArrayOutputFile'</code>
Za zadani genom u FASTA formatu unutar ulazne datoteke se izgradi sufiksno polje i pohrani u zadanu datoteku.

<code>mapper map 'fastaFile' 'suffixArrayOutputFile' 'reads' 'resultOutputFile'</code>
Na zadani referentni genom predan u FASTA formatu, korištenjem prethodno izgrađenog sufiksnog polja se mapiraju kratka očitanja pohranjena u FASTQ formatu unutar predane datoteke. Rezultat u SAM formatu se zapisuje u zadanu izlaznu datoteku.

 
<h4>Dodatno o LCSk matrici:<h4>
1. G. Benson, A. Levy, and B. Shalom, “Longest common subsequence in k length
   substrings,” <br>
2. Filip Pavetić, Goran Žužić, Mile Šikić "A simplified algorithm for the efficient computation of LCSk"
