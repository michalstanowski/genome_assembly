## Opis ogólny
Ten projekt implementuje assembler genomowy oparty o **grafy de Bruijna (DBG)**. Pipeline zawiera:
- **iteracyjną korekcję błędów** na podstawie częstości k‑merów,
- **czyszczenie grafu** (usuwanie ślepych odnóg w grafie, usuwanie bąbli, przycinanie słabych gałęzi),
- **ekstrakcję unitigów**,
- **wydłużanie kontigów poprzez "read-threading"**,
- **scalanie kontigów**,
- **usuwanie redundantnych i niejednoznacznych kontigów**,
- **polishing** metodą dopasowań opartych o seedy

Wejście: FASTA z readami  
Wyjście: FASTA z kontigami

Uruchomienie:
```bash
./assembly input_reads.fasta output_contigs.fasta
```

---

## Podejście krok po kroku

### 1) Wczytanie readów i zapis FASTA
- `load_reads_fasta()` wczytuje sekwencje DNA
- `write_fasta()` zapisuje kontigi w formacie FASTA

### 2) Histogram k‑merów i iteracyjna korekcja błędów
Buduję histogram k‑merów (`kmer_hist`) i lokalnie poprawiam „słabe” k‑mery w odczytach:
- k‑mer jest *słaby* jeśli jego liczność `<= weak`,
- dla słabego k‑meru testuję sąsiadów w odległości Hamminga 1 (`neighbors1mm`),
- opcjonalnie (dla bardzo zaszumionych danych) dopuszczam 2‑mm (`neighbors2mm`),
- podmieniam k‑mer tylko wtedy, gdy najlepszy kandydat jest wyraźnie lepiej wsparty:
  - poprawa co najmniej o `min_improve`,
  - przechodzi test `ratio` względem orgyinalego k-meru,

Funkcja `iterative_correction()` wykonuje te czynności iteracyjnie.

### 3) Wykrywanie „noisy” danych
`estimate_noisy()` ocenia szum na podstawie udziału singletonów 21‑merów.
Jeśli dużo 21‑merów ma count=1, dane uznaję za zaszumione i stosuję agresywniejsze korekcję i czyszczenie.

### 4) Stosowanie wielu różnych wartości k i ogólny zarys pipeline'u
Assembler działa w kilku krokach dla rosnących k (np. 21 → 27 → 31 → 35 → 41 dla danych niezaszumionych). Kroki wymienione poniżej zostaną lepiej wyjaśnione w następnych punktach README.
Dla każdego `(k, min_len)`:
1. ewentualna kolejna korekcja dla danego k (głównie dla małych k),
2. budowa histogramu k‑merów na aktualnym zbiorze sekwencji,
3. wybór progu "solid" (`choose_solid_threshold`), czyli progu określającego czy k-mer nie jest błędem,
4. budowa DBG tylko z "solid" k‑merów,
5. czyszczenie grafu,
6. ekstrakcja unitigów,
7. wydłużanie kontigów,
8. scalenia kontigów,
9. redukcja redundancji i niejednoznaczności,
10. wynik z tego kroku staje się wejściem do następnego kroku k.

Intuicja, którą kierowano się tworząc program:
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

### 6) Upraszczanie grafu (graph cleaning)
Stosowane heurystyki:

**a) `prune_weak_branches`**  
Gdy w węźle jest rozgałęzienie, usuwam wyjścia wyraźnie słabsze od najsilniejszego (korzystając z coverage ratio).

**b) `remove_tips`**  
Usuwam krótkie „tipy” – ślepe odnogi zakończone dead-endem, jeśli ich średnie pokrycie jest niskie.

**c) `remove_bubbles`**  
Wykrywam proste „bąble” (dwie alternatywne ścieżki od tego samego startu do tego samego końca) i usuwam słabszą, jeśli różnica pokrycia jest wystarczająca.

To ogranicza błędy, artefakty i alternatywne ścieżki powodowane szumem.

### 7) Ekstrakcja unitigów
`extract_unitigs()` przechodzi po grafie i składa maksymalne odcinki niejednoznaczności:
- startuje w węzłach, które nie są „1-in-1-out”,
- idzie do przodu dopóki węzły pośrednie są 1-in-1-out,
- zamienia ścieżkę (k‑1)‑merów na sekwencję DNA,
- zostawia tylko sekwencje o długości >= `min_len`.

### 8) Wydłużanie kontigów przez read-threading
Buduję mapy:
- `follow_map[(k-1)-mer][base]` — jak często dana baza występuje **po** (k‑1)‑merze,
- `prev_map[(k-1)-mer][base]` — jak często dana baza występuje **przed** (k‑1)‑merem.

`extend_contig()` próbuje wydłużyć kontig:
- w prawo: na podstawie ostatnich (k‑1) znaków wybiera „pewną” bazę (na podstawie mapy),
- w lewo: analogicznie dla prefiksu,
- zatrzymuje się, gdy wybór staje się niepewny albo k‑mer przestaje być solid.

### 9) Scalanie kontigów przez unikalne nakładki
`merge_by_unique_overlaps(contigs, ovl=k-1)` scala kontigi, jeśli:
- sufiks jednego kontigu pasuje do prefiksu dokładnie jednego innego kontigu (unikalność),
- dzięki temu unikam niejednoznacznych zlepień.

### 10) Redukcja redundancji i kontrola wieloznaczności
Aby nie wypisywać dziesiątek bardzo podobnych kontigów:
- `dedup_contained()` usuwa duplikaty i kontigi zawarte w dłuższych,
- `cap_similar_prefix/suffix/middle()` ogranicza liczbę kontigów o podobnych fragmentach,
- `repetitiveness_ok()` odrzuca bardzo powtarzalne kontigi (mało unikalnych 31‑merów),
- `select_nonredundant_by_kmers()` wybiera podzbiór kontigów, które wnoszą nowe, wspierane k‑mery.

### 11) Polishing (korekta konsensusowa kontigów)
`polish_contigs()` wykonuje kilka rund:
- buduje indeks seedów (`build_seed_index`, np. seed=19) z odczytów,
- dla każdego seeda w kontigu próbuje dopasować odczyt,
- filtruje dopasowania na podstawie mismatchy,
- zbiera głosy per pozycja i poprawia nukleotyd, jeśli ma większość i wystarczająco głosów.

---
