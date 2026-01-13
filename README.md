README

## Opis ogólny
Ten projekt implementuje lekki assembler genomowy oparty o **wielokrotne k (multi-k) grafy de Bruijna (DBG)**. Pipeline zawiera:
- **iteracyjną korekcję błędów** na podstawie częstości k‑merów,
- **czyszczenie grafu** (usuwanie tipów, pękanie bąbli/bubbles, przycinanie słabych gałęzi),
- **ekstrakcję unitigów**,
- **wydłużanie kontigów “read-threading”** (głosowanie baz na końcach na podstawie histogramu),
- **scalanie kontigów przez unikalne nakładki**,
- **filtry anty‑redundancyjne / anty‑wieloznaczności**,
- **polishing** (poprawki konsensusowe) metodą dopasowań opartych o seedy.

Wejście: FASTA z readami  
Wyjście: FASTA z kontigami

Uruchomienie:
```bash
python3 assembly.py input_reads.fasta output_contigs.fasta
```

---

## Podejście krok po kroku (high level)

### 1) Wczytanie readów i zapis FASTA
- `load_reads_fasta()` wczytuje sekwencje DNA (nagłówki `>` są pomijane), wymusza wielkie litery.
- `write_fasta()` zapisuje kontigi w formacie FASTA i zawija linie do 80 znaków.

### 2) Histogram k‑merów i iteracyjna korekcja błędów
Budujemy histogram k‑merów (`kmer_hist`) i lokalnie poprawiamy „słabe” k‑mery w readach:
- k‑mer jest *słaby* jeśli jego liczność `<= weak`,
- dla słabego k‑meru testujemy sąsiadów w odległości Hamminga 1 (`neighbors1mm`),
- opcjonalnie (dla bardzo szumowych danych) dopuszczamy 2‑mm (`neighbors2mm`),
- podmieniamy k‑mer tylko wtedy, gdy najlepszy kandydat jest wyraźnie lepiej wsparty:
  - poprawa co najmniej o `min_improve`,
  - przechodzi test `ratio` względem oryginału,
  - jest lepszy od drugiego najlepszego sąsiada.

Funkcja `iterative_correction()` wykonuje to kilka razy (iteracje), bo po pierwszej rundzie histogram robi się „czystszy”, a kolejne rundy są skuteczniejsze.

### 3) Wykrywanie „noisy” danych
`estimate_noisy()` ocenia szum na podstawie udziału singletonów 21‑merów.
Jeśli dużo 21‑merów ma count=1, dane uznajemy za szumowe i stosujemy agresywniejszą korekcję / czyszczenie.

### 4) Harmonogram multi‑k (kilka wartości k)
Assembler działa w kilku krokach dla rosnących k (np. 21→27→31→35→41 albo krótszy harmonogram dla noisy):
Dla każdego `(k, min_len)`:
1. ewentualnie drobna korekcja dla danego k (głównie dla małych k),
2. histogram k‑merów na aktualnym zbiorze sekwencji,
3. wybór progu „solid” (`choose_solid_threshold`),
4. budowa DBG tylko z solid k‑merów,
5. czyszczenie grafu,
6. ekstrakcja unitigów,
7. wydłużanie kontigów (read-thread),
8. scalenia przez nakładki,
9. filtry i redukcja redundancji,
10. wynik z tego kroku staje się wejściem do następnego kroku k.

Intuicja:
- małe k pomaga połączyć graf mimo błędów i niskiego pokrycia,
- większe k lepiej rozwiązuje powtórzenia i zwiększa jednoznaczność.

### 5) Budowa grafu de Bruijna z solid k‑merów
`build_dbg_from_kmers(kmer_counts, solid)`:
- węzły to (k‑1)‑mery (reprezentowane jako stringi),
- krawędź `u -> v` odpowiada k‑merowi: `u = km[:-1]`, `v = km[1:]`,
- wagi krawędzi to liczności k‑merów (pokrycie),
- struktury:
  - `outC[u][v] = w` (wyjścia),
  - `inC[v][u] = w` (wejścia).

Uwaga: funkcja nie potrzebuje jawnie parametru `k`, bo długość (k‑1) wynika z długości klucza `km`.

### 6) Upraszczanie grafu (graph cleaning)
Stosowane heurystyki:

**a) `prune_weak_branches`**  
Gdy w węźle jest rozgałęzienie, usuwamy wyjścia wyraźnie słabsze od najsilniejszego (coverage ratio).

**b) `remove_tips`**  
Usuwamy krótkie „tipy” – ślepe odnogi zakończone dead-endem, jeśli ich średnie pokrycie jest niskie.

**c) `remove_bubbles`**  
Wykrywamy proste „bąble” (dwie alternatywne ścieżki od tego samego startu do tego samego końca) i usuwamy słabszą, jeśli różnica pokrycia jest wystarczająca.

To ogranicza błędy, artefakty i alternatywne ścieżki powodowane szumem.

### 7) Ekstrakcja unitigów
`extract_unitigs()` przechodzi po grafie i składa maksymalne odcinki niejednoznaczności:
- startujemy w węzłach, które nie są „1-in-1-out”,
- idziemy do przodu dopóki węzły pośrednie są 1-in-1-out,
- zamieniamy ścieżkę (k‑1)‑merów na sekwencję DNA,
- zostawiamy tylko sekwencje o długości >= `min_len`.

### 8) Wydłużanie kontigów przez read-threading (głosowanie baz)
Budujemy mapy:
- `follow_map[(k-1)-mer][base]` — jak często dana baza występuje **po** (k‑1)‑merze,
- `prev_map[(k-1)-mer][base]` — jak często dana baza występuje **przed** (k‑1)‑merem.

`extend_contig()` próbuje wydłużyć kontig:
- w prawo: na podstawie ostatnich (k‑1) znaków wybiera „pewną” bazę (dominujący głos),
- w lewo: analogicznie dla prefiksu,
- zatrzymuje się, gdy wybór staje się niepewny albo k‑mer przestaje być solid.

### 9) Scalanie kontigów przez unikalne nakładki
`merge_by_unique_overlaps(contigs, ovl=k-1)` scala kontigi, jeśli:
- sufiks jednego kontigu pasuje do prefiksu dokładnie jednego innego kontigu (unikalność),
- dzięki temu unikamy niejednoznacznych zlepień.

### 10) Redukcja redundancji i kontrola wieloznaczności
Aby nie wypisywać dziesiątek bardzo podobnych kontigów:
- `dedup_contained()` usuwa duplikaty i kontigi zawarte w dłuższych,
- `cap_similar_prefix/suffix/middle()` ogranicza liczbę kontigów o podobnych fragmentach (buckety),
- `repetitiveness_ok()` odrzuca bardzo powtarzalne kontigi (mało unikalnych 31‑merów),
- `select_nonredundant_by_kmers()` wybiera podzbiór kontigów, które wnoszą nowe, wspierane k‑mery.

To zwykle pomaga w ewaluacji, bo ogranicza multi‑mapping i „szumowe” wyjście.

### 11) Polishing (korekta konsensusowa kontigów)
`polish_contigs()` wykonuje kilka rund:
- buduje indeks seedów (`build_seed_index`, np. seed=19) z readów,
- dla każdego seeda w kontigu próbuje dopasować read (offset),
- filtruje dopasowania zbyt wieloma mismatchami,
- zbiera głosy per pozycja i poprawia bazę, jeśli ma większość i wystarczająco głosów.

Polishing poprawia identyczność (identity), zwłaszcza gdy wcześniejsze etapy zostawiły pojedyncze błędy.

---

## Dlaczego czasem „ucina” liczbę kontigów?
W tej wersji kodu **nie ma twardego `[:80]` na końcu** (to było w jednej z wcześniejszych wersji).
Gdyby limit był włączony, służyłby głównie do:
- ograniczenia redundancji,
- zmniejszenia ryzyka multi‑mappingu,
- skrócenia czasu ewaluacji i rozmiaru wyjścia.

Obecnie kontigi są tylko sortowane po długości i wszystkie są zapisywane (po filtrach).

---

## Ograniczenia / uwagi
- Pipeline nie uwzględnia reverse-complement (jedna orientacja).
- Bubble popping jest „prosty” (nie obsłuży złożonych struktur).
- Polishing to seed-based pileup, nie pełne globalne wyrównanie.
- Progi/parametry są heurystyczne i mogą wymagać strojenia pod dane.
