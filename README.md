projekt-mapper
==============

// POTREBNO AZURIRATI OPIS

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

<p>
<code>mapper -m index 'fastaFile' 'suffixArrayOutputFile'</code><br>
Za zadanu sekvencu u FASTA formatu unutar ulazne datoteke se izgradi sufiksno polje i pohrani u zadanu datoteku.
</p>
<p>
<code>mapper -m map [options] 'fastaFile' 'suffixArrayOutputFile' 'readsInputFASTQ' 'resultOutputFile'</code><br>
Na zadani referentni genom predan u FASTA formatu, korištenjem prethodno izgrađenog sufiksnog polja se mapiraju kratka očitanja pohranjena u FASTQ formatu unutar predane datoteke. Rezultat u SAM formatu se zapisuje u zadanu izlaznu datoteku.

Options:
<ul>
<li>-t N      N is thread number. [default: number of cores]</li>
<li>-k N      N is seed length. [default: 15]</li>
<li>-l lowerLimit     lowerLimit is minimum required number of hits. [default: 45 ]</li>
<li>-h upperLimit    upperLimit is maximum allowed number of hits. [default: 20 ]</li>
<li>-kf factor     factor float value that represents minimum ratio bestScore/score of reads that are being kept. [default: 1.2 ]</li>
<li>--pos N     N is maximum allowed number of positions that are kept per read [default: 80 ]</li>
</ul>

Example: 
<p>
Create index
<code>./mapper -m index caenorhabditis_elegans.fa caenorhabditis_elegans.fa.suffa </code>
</p>

<p>
Map reads
<code>./mapper -m map -t 4 -k 10 -l 20 -h 45 caenorhabditis_elegans.fa caenorhabditis_elegans.fa.suffa caenorhabditis_elegans/reads.fa MAPPER-10-45-20-v1.sam</code>
</p>



 
<h4>Dodatno o LCSk matrici:</h4>
1. G. Benson, A. Levy, and B. Shalom, “Longest common subsequence in k length
   substrings,” <br>
2. Filip Pavetić, Goran Žužić, Mile Šikić "A simplified algorithm for the efficient computation of LCSk"
