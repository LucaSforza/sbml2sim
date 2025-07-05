#set page(header: align(right)[Luca Sforza 2050030], numbering: "1")

#set heading(numbering: "1.")

#set par(
  first-line-indent: 1em,
  spacing: 0.65em,
  justify: true,
)

#set figure(supplement: "Plot", numbering: none)

#set ref(supplement: none)

#set math.equation(numbering: "(1)")


#align(center, text(17pt)[
  *Simulazione del ripiegamento proteico mediato da chaperoni nel cancro* 
])

#image("R-HSA-392499.svg")

= Introduzione

La simulazione quantitativa di processi biologici complessi richiede modelli dinamici basati su leggi cinetiche con parametri precisi, come le costanti cinetiche delle reazioni.

Tuttavia, molti modelli biologici disponibili, come quelli di Reactome @milacic2024reactome,
sono qualitativi e non forniscono questi dati essenziali, limitando la possibilità di simulazioni realistiche.

In questo lavoro, ci si concentra sulla trasformazione di un modello qualitativo di ripiegamento proteico mediato da chaperoni @reactome_protein_folding
in un modello quantitativo utilizzabile per simulazioni dinamiche.

Per superare la mancanza di dati sperimentali sulle costanti cinetiche, si utilizza un approccio di ottimizzazione black box per stimare questi parametri, vincolando il modello a dati di concentrazione proteica tipici del cancro al seno.

= Simulazione del modello

== Formato SBML Livello 3 Versione 1 e descrizione del file SBML relativo al pathway di riferimento

I modelli di Reactome sono scaricabili con il formato SBML (System Biology Markup Language) basato su XML @hucka2003sbml.

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

Invece il moto è descritto dalle leggi cinetiche delle reazioni, che a loro volta descrivono la velocità di una reazione @klipp2009systems[cap. 3].

Sia $R = {1,...,n}$ l'insieme degli identificatori delle reazioni, sia $S$ la concentrazione di una singola specie. $v_i$ è la velocità della reazione $i in R$. $R_S subset.eq R$ è l'insieme delle reazioni dove $S$ è reagente. $P_S subset.eq R$ è l'insieme delle reazioni dove $S$ è prodotto. $n^S_i$ è la stechiometria di $S$ nella reazione $i in R$. 

Allora:

$
  (d S)/(d t) = sum_(i in P_S) n^S_i v_i - sum_(j in R_S) n^S_j v_j
$

== Leggi Cinetiche

Dato che i modelli di Reactome sono qualitativi e non quantitativi mancano totalmente le leggi cinetiche. Le legge utilizzata è la mass action rule, la velocità di una reazione è proporzionale alla concentrazione dei reagenti. Ad essa è stata aggiunta la _hill function_ per modellare i modificatori delle reazioni.
Quindi la velocità $v$ di una reazione $R$ è definita come segue:

Siano $A_1, A_2, ..., A_n$ le specie reagenti, $n_i$ la stechiometria del reagente $i$ e $M_1, M_2, ..., M_m$ le specie modificatrici.

$
  v = product_(i=1)^m H(M_i) dot K_R dot product_(i=1)^n A_i^(n_i)
$
#pagebreak()

$K_R$ rappresenta la costante cinetica specifica della reazione $R$ ed è aggiunta come parametro del modello.

$H(S)$ è la hill function che è definita come segue:

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

$accent(x, tilde)(t, theta)$ = $x(t, theta) + epsilon$ dove $epsilon$ è un disturbo casuale. 
Il disturbo avviene perché alcune specie in input non si conoscono le concentrazioni medie,
quindi viene assegnato un valore casuale in un certo range realistico.

Quindi le costanti cincetiche da trovare devono minimizzare la loss function per un qualsiasi valore
per le specie in input di cui non si hanno i dati.

Inoltre le costanti cinetiche, per essere realistiche, devono portare lo stato del sistema in uno stato di equilibrio.

Dobbiamo quindi definire una funzione di utilità che penalizzi i sistemi che non raggiungono uno stato stazionario.

Sia $accent(m, tilde)_i (t, theta)$ la concentrazione media della specie $i$ al tempo $t$ con i parametri $theta$.

Introduciamo un iper-parametro $phi in [0,1]$ che rappresenta la frazione dell’orizzonte temporale $T$ considerata per valutare la stabilità.

Definiamo quindi:

$
  cal(L)_1(theta) =  sum_(i=1)^n (accent(m, tilde)_i (phi dot T, theta) - accent(m, tilde)_i (T, theta))^2
$

Per i test, ho scelto $phi approx 0.80$, in modo da valutare la variazione delle concentrazioni medie tra la fase finale e quella immediatamente precedente della simulazione. Questo permette di ignorare le fluttuazioni iniziali dovute alle condizioni iniziali e concentrarsi sulla stabilità asintotica del sistema.

Adesso dobbiamo vincolare il valore delle specie di cui si conosce la concentrazione media.

Sia: $
cal(D) = {(S_i, y) | S_i "si conoscono le concentrazioni" and "y è la concentrazione mol/L"}$

$
  cal(L)_2(theta) = sum_((S_i, y) in cal(D)) (log_10 (accent(x,tilde)_i (T, theta)) - log_10 (y))^2
$

Ho adottato un errore quadratico sulla scala logaritmica per la funzione di loss $LL_2$ al fine di confrontare in modo coerente le concentrazioni simulate con i dati osservati. Questo approccio consente di trattare in maniera naturale quantità che variano su piu' ordini di grandezza.

Durante la simulazione possono verificarsi errori di integrazione numerica, tipicamente quando le costanti cinetiche assumono valori troppo elevati, causando una variazione troppo rapida delle concentrazioni delle specie che tendono rapidamente a $-infinity$ o $+infinity$.

Le costanti cinetiche scelte non devono avere questa caratteristica quindi la terza funzione di loss sarà o 0 oppure $+infinity$.

$
  cal(L)_3(theta) = cases(
    +infinity "se il sistema ha ottenuto errori di integrazione numerica",
    0 "altrimenti"
  )
$

La loss function finale è:

$
cal(L)(theta) = S dot [ space p dot cal(L)_1(theta) + (1 - p) dot cal(L)_2(theta) space] + cal(L)_3(theta)
$

Dove $p in [0,1]$ è un iper-parametro che bilancia l'importanza relativa tra la stabilità del sistema ($cal(L)_1$) e l'aderenza ai dati sperimentali ($cal(L)_2$). Invece $S in RR$ è il fattore di scala.

Valori di $p$ prossimi a 1 privilegiano la stabilità, mentre valori vicini a 0 danno maggiore importanza alla corrispondenza con i dati sperimentali. La funzione $cal(L)_3$ agisce come vincolo rigido, penalizzando con $+infinity$ le soluzioni che generano errori numerici, e pertanto non necessita di un coefficiente di ponderazione.

Per questo progetto ho scelto $S = 10^4$ e $p = 0.1$, quindi ho privileggiato l'aderenza ai dati sperimentali.

#import "@preview/lovelace:0.3.0": *

== Funzionamento dell'ottimizzazione Black-Box di Nevergrad

La suite di ottimizzatori usati per questo progetto è Nevergrad @nevergrad, sviluppato da Meta.

Nevergrad utilizza una vasta gamma di ottimizzatori utilizzabili. Uno di questo è *Wizard* che opera come meta-euristica su quale algoritmo di ottimizzazione va usato, selezionando anche gli iper-parametri.

Inoltre gli ottimizzatori di Nevegrad si basano sul patter _ask and tell_.

Ossia ogni ottimizzatore ha un intefaccia in cui permettono di richiedere dei parametri
che rappresenta il tentativo di ottimizzare la funzione. ($theta <- "optimizer"."ask()"$).

Invece la _tell_ permette di informare l'ottimizzatore la _loss_ dei parametri scelti.($"optimizer"."tell"(theta, cal(L)(theta))$)

Per valutare le performance lo si può fare con $LL("optimizer"."recommend")$.

Per ottimizzare la funzione è stato usato il seguente algoritmo.
#figure(
pseudocode-list[
  + *for* $i$ *in* range(*budget*) *do*
    + $theta <- "optimizer"."ask"()$;
    + $"optimizer"."tell"(theta, cal(L)(theta))$
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
    + Values $<- {(theta, cal(L)(theta)) | theta in "DecisionSet"}$ $"#"$ calcolato in parallelo
    + $"optimizer"."tell"("Values")$
  + *end for*;
]
)

#pagebreak()

Questo algoritmo ha due vantaggi rispetto al precedente:

+ La funzione di _loss_ può essere parallelizzata. RoadRunner (simulatore di sistemi biologici) permette di eseguire in parallelo piu' modelli.
+ *Wizard* può scegliere un algoritmo diverso per ottimizzare che sfrutta il fatto di poter richiedere più parametri contemporaneamente rispetto dal grado di parallelismo (numero di soluzione che possono essere calcolate contemporaneamente).

== Alcuni algoritmi utilizzati da Nevergrad

Tra gli algoritmi che *Wizard* può scegliere ce ne sono veramente tanti, ma ecco una breve descrizione dei principali algoritmi.

=== Random Search

La random search @wikipedia-randomsearch è uno degli algoritmi di ottimizzazione black box più semplici. All'inizio seleziona un parametro casuale nello spazio delle soluzioni. Successivamente, ad ogni iterazione, genera nuovi parametri casuali all'interno di un'ipersfera di raggio $r$ centrata sull'attuale soluzione migliore (raccomandation). Se la loss calcolata sui nuovi parametri è inferiore a quella corrente, la raccomandation viene aggiornata con i nuovi valori trovati.

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


=== OnePlusOne 

L'algoritmo OnePlusOne Evolution Strategy @one_plus_one_es è una strategia euristica di ottimizzazione iterativa basata su mutazioni casuali.

L'idea è molto semplice: si parte da una soluzione casuale iniziale, e ad ogni iterazione viene generata una nuova soluzione mutata a partire da quella attuale (chiamato anche _parent_). Se la nuova soluzione ha una loss inferiore rispetto a quella precedente, allora viene accettata come nuova raccomandation.

La mutazione è un disturbo sulle variabili di decisione $cal(N) (0,sigma)$ e $sigma$ cambia dinamicamente durante l'esecuzione dell'algoritmo in base al successo delle mutazioni precedenti, secondo una regola euristica chiamata 1/5th success rule.

Secondo questa regola, se più di 1 mutazione su 5 viene accettata (cioè migliora la loss), allora l'algoritmo aumenta il passo di mutazione $sigma$ per esplorare più velocemente. Se invece meno di 1 su 5 migliora, $sigma$ viene ridotto, per favorire l'esplorazione locale.

#figure(
pseudocode-list[
  + ask(*self*: _optimizer_) $->$ _Parameter_:
    + *if* *self*.raccomandation *is* *None* *then* *return* parametro casuale
    + *return* parametro mutato da *self*.raccomandation con deviazione standard $sigma$
  + tell(*self*: _optimizer_, $theta$: _Parameter_, loss: _Real_):
    + *if* loss < *self*.raccomandation.loss *then*
      + *self*.raccomandation = $theta$
      + aggiungi 1 alla finestra di successi
    + *else*:
      + aggiungi 0 alla finestra di successi
    + *end if*;
    + *if* è il momento di aggiornare $sigma$ *then*
      + *if* frequenza di successi > 1/5 *then* aumenta $sigma$
      + *else* diminuisci $sigma$
    + *end if*;
]
)
)

L'algoritmo OnePlusOne è particolarmente utile quando lo spazio dei parametri è continuo e di dimensione moderata, ed è stato ampiamente usato nel contesto dell'ottimizzazione evolutiva. Il suo vantaggio principale è la capacità di adattarsi automaticamente alla scala del problema, migliorando progressivamente l'efficienza della ricerca.


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

#pagebreak()

#bibliography("refs.bib", title: "Bibliografia")
