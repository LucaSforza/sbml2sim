#set page(header: align(right)[Luca Sforza 2050030], numbering: "1")

#set heading(numbering: "1.")

#set par(
  first-line-indent: 1em,
  spacing: 0.65em,
  justify: true,
)

#set figure(supplement: "Plot")

#set ref(supplement: none)

#set math.equation(numbering: "(1)")

#align(center, text(17pt)[
  *Simulazione del ripiegamento proteico mediato da chaperoni nel cancro* 
])

#image("R-HSA-392499.svg")

= Introduzione

La simulazione quantitativa di processi biologici complessi richiede modelli dinamici basati su leggi cinetiche con parametri precisi, come le costanti cinetiche delle reazioni.

Tuttavia, molti modelli biologici disponibili, come quelli di Reactome, // TODO: riferimento reactome
sono qualitativi e non forniscono questi dati essenziali, limitando la possibilità di simulazioni realistiche.

In questo lavoro, ci si concentra sulla trasformazione di un modello qualitativo di ripiegamento proteico mediato da chaperoni // TODO: riferimento
in un modello quantitativo utilizzabile per simulazioni dinamiche.

Per superare la mancanza di dati sperimentali sulle costanti cinetiche, si utilizza un approccio di ottimizzazione black box per stimare questi parametri, vincolando il modello a dati di concentrazione proteica tipici del cancro al seno.

= Simulazione del modello

== Formato SBML Livello 3 Versione 1 e descrizione del file SBML relativo al pathway di riferimento

I modelli di Reactome sono scaricabili con il formato SBML (System Biology Markup Language) basato su XML. // TODO: riferimento a SBML

Per questo progetto è stato utilizzato il livello 3 versione 1, poiché è quello che offre Reactome tra i suoi modelli.

La parte più esterna di un file SBML è rappresentata dall'oggetto SBML:

```XML
<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">
  ...
  <model ...> ...
  </model>
</sbml>
```

Esso ha come unica funzione quello di specificare il livello e la versione, ad esso è associato un ogetto di tipo "model"
in cui è contenuto la descrizione del modello biologico.

Un modello è semplicemente un container per altri oggetti, nel pathway di riferimento di questo progetto
sono presenti solo tre tiplogie di oggetti: I compartienti, Le specie e le reazioni.

Un compartimento in SBML rappresenta uno spazio delimitato dove sono localizzate le specie.

Anche se i compartimenti sono opzionali, una specie ha l'obbligo di specificare il compartimento di appartenenza, quindi sono di fatto obbligatori.

I compartimenti  nel pathway di riferimento sono: membrana plasmatica e citosol.
Sono state utilizzate come grandezze quantitative per questi compartimenti quelle del cancro al seno.

Una specie in SBML rappresenta un insieme di entità indistinguibili tra di loro, possono partecipare a reazioni e sono localizzati in uno specifico compartimento.

Una reazione in SBML descrive ogni tipo di processo che cambia la quantità di una o più specie. Una reazione in SBML necessariamente deve definire le sue proprietà strutturali, ovvero specificare i reagenti e/o i prodotti (volendo anche i modificatori). Una reazione può (ma non è obbligata) ad avere pure una sua descrizione quantitativa della reazione, ovvero una legge cinetica.

Le leggi cinetiche descrivono la velocità della reazione e da quelle si può ricavare il moto del sistema.



== Differenza tra Modelli Qualitativi e Quantitativi <qual>

Reactome offre modelli biologici qualitativi e non quantitativi. Questo vuol dire che per quanto riguarda le reazioni mancano le leggi cinetiche, per i vari compartimenti non è specificato il loro volume, non sono specificate le unità di misura e non sono specificate le concentrazioni iniziali delle specie.

Un modello qualitativo ha come obiettivo non quello di essere simulabile, ma quello di descrivere il modello biologico.

SBML può essere usato anche per descrivere modelli quantitativi.

Un modello biologico può essere visto come un sistema dinamico non lineare tempo invariante.

Lo stato del sistema se $n$ sono le specie è $x = mat(S_1;...;S_n) in RR^n$.

Dove $S_i$ è la concentrazione della specie $i$.

Invece il moto è descritto dalle leggi cinetiche delle reazioni, che a loro volta descrivono la velocità di una reazione.

Sia $S$ la concentrazione di una singola specie. $v_1,...,v_k$ è la velocità delle reazioni dove $S$ è reagente e $v'_1,...,v'_m$ è la velocità delle reazioni dove $S$ è prodotto. // TODO: riferimento system biology

Allora: // TODO: riferimento system biology

$
  (d S)/(d t) = sum_(i=1)^m v'_i - sum_(j=1)^k v_i
$

== Leggi Cinetiche

Dato che i modelli di Reactome sono qualitativi e non quantitativi mancano totalmente le leggi cinetiche. Le legge utilizzata è quella di Michelis-Mentent. // TODO: non è la mass action?
Quindi la velocità $v$ di una reazione $R$ è definita come segue:

Siano $A_1, A_2, ..., A_n$ le specie reagenti, $n_i$ la stechiometria del reagente $i$ e $M_1, M_2, ..., M_m$ le specie modificatrici.
// TODO: riferimento system biology mass action
$
  v = product_(i=1)^m H(M_i) dot K_R dot product_(i=1)^n A_i^(n_i)
$

$K_R$ rappresenta la costante cinetica specifica della reazione $R$ ed è aggiunta come parametro del modello.

$H(S)$ è la hill function che è definita come segue: // TODO: riferimento hill function

$
  H( S) = cases(
    (S^h)/(K_(a,R)^h + S^h) "se "S" è attivatore",
    (K_(i,R)^h)/(K_i^h + S^h) "se "S" è un inibitore" ,
  )
$

Tramite i file SBML posso ottenere come informazione se $S$ è un inibitore oppure no dall'identificatore SBO (System Biology Ontology).

$K_(a,R)$ e $K_(i,R)$ sono delle nuove costanti aggiunte ai parametri del modello.

== Calcolo medie del valore delle concentrazioni

Per stimare le costanti cinetiche è necessario sapere la concentrazione media delle specie alla fine della simulazione, questo dato serve per due motivi:

+ Verificare che la concentrazione media è quella che ci si aspettava rispetto ai dati sperimentali.
+ Verificare che il sistema sia stabile.

Per questo motivo, ho aggiunto un modulo che introduce nel modello SBML dei parametri variabili che rappresentano la concentrazione media di ciascuna specie.

Sia $x_S$ la concentrazione media della specie $S$. Allora la sua dinamica è definita come segue:

$
  (d x_S)/(d t) = (S(t) -x_S (t))/(t + epsilon)
$ <din>

Dove $epsilon$ è un numero sufficientemente piccolo, per gli scopi di questo progetto ho scelto come valore: $10^(-6)$.

L'equazione differenziale @din se simulata con un orizzonte abbastanza grande il valore di $x_S$ combacerà con il vero valore medio della concentrazione $S$.

In più abbiamo che $x_S (0) = S(0)$.

== Input e Output

Un modello biologico normalmente ha delle specie di Input e delle specie di Output. Riconoscerle è abbastanza semplice; Se una specie 
appare solo come reagente delle reazioni allora è un Input. Se una specie appare solo come prodotto allora è un Output.

Pertanto, per ogni specie di input la concentrazione viene mantenuta costante durante la simulazione, mentre per ogni specie di output la concentrazione iniziale è posta a zero e rimane costante, simulando così la rimozione continua dell'output dal sistema.

= Stima delle costanti cinetiche

Per le costanti cinetiche non abbiamo dati sperimentali su cui affidarci, quindi vanno stimate.

Un modo per farlo è definire una *loss function* e minimizzarla ottenendo cosi' dei parametri realistici per il sistema.

Siano $theta$ le costanti cincetiche del sistema (ovvero i paramentri).

Dovremmo dare un dominio a queste variabili, anche perché l'ottimizzatore riesce ad essere piu' efficiente se i parametri sono *bounded*.

Se $theta_i$ è una singola costante cinetica, allora un range di valori ragionevole è $theta_i in [10^(-6), 10^6]$.

Tuttavia, questo approccio presenta una problema: le costanti cinetiche possono variare su diversi ordini di grandezza. Di conseguenza, un intervallo ampio come $[10^(-6), 10^6]$ può causare problemi all'ottimizzatore, il quale tende a esplorare maggiormente le regioni dell'intervallo con valori elevati, trascurando invece le zone vicine allo zero. Questo squilibrio nella distribuzione dei punti esplorati può compromettere l'efficacia della ricerca.

Per superare questa limitazione, è possibile adottare un'ottimizzazione in scala logaritmica.

Invece di ottimizzare direttamente i valori delle costanti cinetiche, si ottimizzano i loro esponenti in base 10.

In questo modo, il problema viene riformulato come la ricerca di un vettore $theta in [-6, 6]^d$, dove
$d$ è il numero di costanti cinetiche e ciascuna costante è poi ricostruita come $k_i = 10^(theta_i)$.

#pagebreak()

== Loss function

L'ottimizzatore deve indovinare quindi un vettore $theta in [-6, 6]^d$ che minimizza una loss function da definire.

La loss function deve codificare l'utilità dei parametri scelti.

I requisiti per le costanti cinetiche da scegliere sono i seguenti:

+ Il sistema deve terminare in uno stato di stabilità.
+ Il valore delle proteine di cui si conoscono le concentrazione medie devono combaciare con il valore delle concentrazioni simulate.

=== Modellazione della Loss function

Sia $n$ il numero di specie nel modello.

Sia $theta in [-6, 6]^d subset RR^d$ il vettore degli esponenti delle costanti cinetiche:

$
  k_i = 10^(theta_i)
$

Dove $k_i$ è la costante cinetica $i$.

Definiamo una funzione di utilità delle costanti cinetiche.

Sia $T$ l'orizzonte della simulazione.

$x(t, theta)$ lo stato del sistema al tempo $t$ con i parametri del sistema $theta$.

$accent(x, hat)(t, theta)$ = $x(T, theta) + epsilon$ dove $epsilon$ è un errore casuale 
L'errore casuale avviene poiché per le specie in input di cui non si conoscono le concentrazioni medie
viene assegnato un valore casuale in un certo range realistico.

Quindi le costanti cincetiche da trovare devono minimizzare la loss function per un qualsiasi valore
per le specie in input di cui non si hanno i dati.
// TODO: dire meglio


Le costanti cinetiche, per essere realistiche, devono portare lo stato del sistema in un punto di equilibrio.

Dobbiamo quindi definire una funzione di utilità che penalizzi i sistemi che non raggiungono uno stato stazionario.

Sia $accent(m, hat)_i (t, theta)$ la concentrazione media della specie $i$ al tempo $t$ con i parametri $theta$.

Introduciamo un iper-parametro $phi in [0,1]$ che rappresenta la frazione dell’orizzonte temporale $T$ considerata per valutare la stabilità.

Definiamo quindi:

$
  LL_1(theta) =  sum_(i=1)^n (accent(m, hat)_i (phi dot T, theta) - accent(m, hat)_i (T, theta))^2
$

Per i test, ho scelto $phi approx 0.80$, in modo da valutare la variazione delle concentrazioni medie tra la fase finale e quella immediatamente precedente della simulazione. Questo permette di ignorare le fluttuazioni iniziali dovute alle condizioni iniziali e concentrarsi sulla stabilità asintotica del sistema.

Adesso dobbiamo vincolare il valore delle specie di cui si conosce la concentrazione media.

Sia: $
DD = {(S_i, y) | S_i "si conoscono le concentrazioni" and "y è la concentrazione mol/L"}$

$
  LL_2(theta) = sum_((S_i, y) in DD) (accent(x,hat)_i (T, theta) - y)^2
$

Durante la simulazione possono verificarsi errori di integrazione numerica, tipicamente quando le costanti cinetiche assumono valori troppo elevati, causando una variazione troppo rapida delle concentrazioni delle specie che tendono rapidamente a $-infinity$ o $+infinity$.

Le costanti cinetiche scelte non devono avere questa caratteristica quindi la terza funzione di loss sarà o 0 oppure $+infinity$.

$
  LL_3(theta) = cases(
    +infinity "se il sistema ha ottenuto errori di integrazione numerica",
    0 "altrimenti"
  )
$

#pagebreak()

La loss function finale è:

$
LL(theta) = p dot LL_1(theta) + (1 - p) dot LL_2(theta) + LL_3(theta)
$

Dove $p$ è un iper-parametro che bilancia l'importanza relativa tra la stabilità del sistema ($LL_1$) e l'aderenza ai dati sperimentali ($LL_2$).

Valori di $p$ prossimi a 1 privilegiano la stabilità, mentre valori vicini a 0 danno maggiore importanza alla corrispondenza con i dati sperimentali. La funzione $LL_3$ agisce come vincolo rigido, penalizzando con $+infinity$ le soluzioni che generano errori numerici, e pertanto non necessita di un coefficiente di ponderazione.

// TODO: dire quale iper-parametro ho scelto per p

#import "@preview/lovelace:0.3.0": *

== Funzionamento dell'ottimizzazione Black-Box di Nevergrad
// TODO: leggere articoli

La suite di ottimizzatori usati per questo progetto è Nevergrad, sviluppato da Meta. // TODO: riferimento articolo

Nevergrad utilizza una vasta gamma di ottimizzatori utilizzabili. Uno di questo è *Wizard* che opera come meta-euristica su quale algoritmo di ottimizzazione va usato, selezionando anche gli iper-parametri.

Inoltre gli ottimizzatori di Nevegrad si basano sul patter _ask and tell_.

Ossia ogni ottimizzatore ha un intefaccia in cui permettono di richiedere dei parametri
che rappresenta il tentativo di ottimizzare la funzione. ($theta <- "optimizer"."ask()"$).

Invece la _tell_ permette di informare l'ottimizzatore la _loss_ dei parametri scelti.($"optimizer"."tell"(theta, LL(theta))$)

Per valutare le performance lo si può fare con $LL("optimizer"."recommend")$.

Per ottimizzare la funzione è stato usato il seguente algoritmo.
#figure(
pseudocode-list[
  + *for* $i$ *in* range(*budget*) *do*
    + $theta <- "optimizer"."ask"()$;
    + $"optimizer"."tell"(theta, LL(theta))$
  + *end for*;
]
)

Ma si potrebbe fare di meglio con il seguente algoritmo (non implementato per questo progetto, solo un idea).

#figure(
pseudocode-list[
  + *for* $i$ *in* range(*budget*) *do*
    + DecisionSet $<- emptyset.rev$
    + *for* $"_"$ *in* range(*parallel degree*) *do*
      + DecisionSet $<-$ DecisionSet $union {"optimizer"."ask"()}$
    + *end for*;
    + Values $<- {(theta, LL(theta)) | theta in "DecisionSet"}$ $"#"$ calcolato in parallelo
    + $"optimizer"."tell"("Values")$
  + *end for*;
]
)

Questo algoritmo ha due vantaggi rispetto al precedente:

+ La funzione di _loss_ può essere parallelizzata. RoadRunner (simulatore di sistemi biologici) permette di eseguire in parallelo piu' modelli.
+ *Wizard* può scegliere un algoritmo diverso per ottimizzare che sfrutta il fatto di poter richiedere più parametri contemporaneamente rispetto dal grado di parallelismo (numero di soluzione che possono essere calcolate contemporaneamente).

#pagebreak()

== Algoritmi utilizzati da Nevergrad

Tra gli algoritmi che *Wizard* può scegliere ce ne sono veramente tanti, ma ecco una breve descrizione dei principali algoritmi e come vengono scelti da Wizard.

=== Random Search // TODO: ref: https://en.wikipedia.org/wiki/Random_search

Questo è il piu' semplice algoritmo di black box utilizzabile. Semplicemente esegue esperimenti sulla funzione
obbiettivo usando parametri casuali e da come raccomandation i parametri che minimizzano la loss function.

#figure(
  pseudocode-list[
    + ask(*self*: _optimizer_) $->$ _Parameter_:
      + *if* *self*.raccomandation *is* *None* *then* *return* parametro casuale
      + Sia *self*.r iper-parametro del modello che è il raggio di una ipersfera
      + *return* parametro casuale all'interno dell'ipersfera di raggio *self*.r con centro *self*.raccomandation
    + tell(*self*: _optimizer_, $theta$: _Parameter_, loss: _Real_):
      + *if* *self*.raccomandation.loss > loss *then*
        + *self*.raccomandation = $theta$
  ]
)


=== Bayesian Optimization // TODO: ref BoTorch

Questa tecnica di ottimizzazione è composta da 2 componenti:
+ Un surrogato probabilistico $f$ rispetto alla vera funzione $f_"true"$ da ottimizzare (generalmente un processo gaussiano)
+ Una funzione di acquisizione $alpha$ che ottimizza il surrogato probabilistico al posto della funzione vera.

// TODO: descrivi i processi gaussiani

Un processo gaussiano è un modello probabilistico utilizzato per approssimare funzioni sconosciute in modo non parametrico.

Sia $DD = {(x_i, y_i) | i in [1, n]}$ il dataset degli $n$ esperimenti osservati, dove $y_i = f_"true" (x_i) + epsilon$ e $epsilon$ rappresenta il rumore.

L'obiettivo è stimare il valore della funzione per un nuovo input $x$, cioè calcolare la distribuzione predittiva $P(y| x, DD)$ che si assume essere una gaussiana.

Un processo gaussiano inoltre assume che ogni insieme finito di punti $(x_1, ..., x_n)$ abbia una distribuzione congiunta gaussiana, specificata da una media e una funzione di covarianza (kernel). Questo permette di stimare la media e la varianza della funzione nei punti non osservati, fornendo sia una previsione che una misura di incertezza.

// TODO: scrivi la ask e la tell

#figure(
  pseudocode-list[
    + ask(*self*: _optimizer_) $->$ _Parameter_:
      + *return* *self*.$alpha$(*self*.GaussianProcess);
    + tell(*self*: _optimizer_, $theta$: _Parameter_, loss: _Real_):
      + *self*.GaussianProcess.$DD$ = *self*.GaussianProcess.$DD union {(theta, "loss")}$
  ]
)


=== OnePlusOne // TODO: ref: https://algorithmafternoon.com/strategies/one_plus_one_evolution_strategy/



#pagebreak()

= Risultati
#figure(
  grid(
    columns: 1,
    rows: 3,
    //gutter: 2mm,
    grid(
      columns: 3,
      image("breast_cancer_cell_species_379537.png"),
      image("breast_cancer_cell_species_379538.png"),
      
    ),
    grid(
      columns: 2,
      image("breast_cancer_cell_species_379539.png"),
      image("breast_cancer_cell_species_379540.png"),
    ),
      image("breast_cancer_cell_species_379546.png", width: 55%)
  ),
   caption: [Concentrazione delle specie vincolate dalla loss function $LL_2$ rispetto ai dati sperimentali],
)

All'inizio della simulazione, i valori delle concentrazioni vengono impostati uguali a quelli ottenuti dai dati sperimentali.
Questo approccio serve a far partire il simulatore già vicino a uno stato di equilibrio.
Partire vicino all'equilibrio rende più semplice e veloce il processo di ottimizzazione, perché il simulatore deve fare meno aggiustamenti per trovare i parametri migliori.

Tuttavia, non tutte le specie chimiche riescono a raggiungere esattamente i valori sperimentali. Questo è normale, perché l'ottimizzazione black box è una tecnica che non assicura l'ottimalità, quindi può non trovare la soluzione perfetta per tutte le specie.

#figure(
  grid(
    columns: 2,
    image("kinetic2.png"),
    image("kinetic_log.png")
  ),
  caption: [Concentrazioni medie di tutte le specie presenti nel sistema: a sinistra scala lineare, a destra scala logaritmica.]
)

Da questo ultimo plot si può visualizzare la stabilità del sistema. All'inizio della simulazione non è stabile, ma arriva velocemente ad un punto di stabilità.

#figure(
  image("log.png", width: 80%),
  caption: [Andamento della funzione di utilità (loss) durante i tentativi dell’ottimizzatore.]
)