#import "@preview/slydst:0.1.4": *
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

#show: slides.with(
  title: "Stima delle costanti cinetiche delle reazioni di un modello biologico",
  // subtitle: "Algoritmo naive per il clustering",
  authors: "Luca Sforza 2050030",
)

#show raw: set block(fill: silver.lighten(65%), width: 100%, inset: 1em)

#set ref(supplement: none)

#set math.equation(numbering: "(1)")

=== Contenuti

#outline()

= Introduzione

== Cos'è un modello?

Un modello è una rappresentazione astratta di un sistema reale, che permette di descriverne il comportamento

#let bent-edge(from, to, ..args) = {
  let midpoint = (from, 50%, to)
  let vertices = (
    from,
    // (from, "|-", midpoint),
    // (midpoint, "-|", to),
    to,
  )
  edge(..vertices, "-|>", ..args)
}

#figure(


  diagram(
  node-stroke: luma(80%),
  edge-corner-radius: none,
  spacing: (30pt, 20pt),

  // Nodes
  node((-1, 0), [*$u$*], name: <a>),
  node((0, 0), [*$h$*], name: <b>),
  node((1, 0), [*$y$*], name: <c>),
  node((0, -1), [$theta$], name: <d>),
  // Edges
  bent-edge(<a>, <b>),
 bent-edge(<b>, <c>),
 bent-edge(<d>, <b>)
)
)

Dati gli input $u$, ad esempio condizioni iniziali o ingressi esterni e i parametri $theta$ , il modello produce come output $y$. In altre parole, il modello descrive la relazione:
$
  y = h(u, theta)
$ // TODO: dire meglio

== Cos'è un modello biologico

Nello specifico un modello biologico è un sistema dinamico di equazioni differenziali ordinarie.

Questo tipo di modello descrive come evolve nel tempo lo stato del sistema in funzione di certi parametri.

Lo stato è rappresentato da un vettore $x(t, theta) in RR^n$, che cambia nel tempo $t$ e dipende da alcuni parametri $theta$.

Ogni elemento $x_i (t, theta)$ rappresenta la concentrazione della specie $i$-esima, cioè un insieme di entitià chimiche considerate indistinguibili tra loro.

Una specie per far parte dello stato deve essere sia reagente che prodotto di almeno una reazione chimica. Se è solo reagente allora è un input, se è solo prodotto allora è un output.

== Reactome

// TODO: da riadattare tutto il documento

Molti modelli biologici sono disponibili, come quelli di *Reactome* @milacic2024reactome database di pathway scritti usando SBML (System Biology Markup Language).

Questi sono modelli sono qualitativi e non descrivono la dinamica del sistema o gli stati iniziali, ma solo la sua descrizione biologica.

Per poter rendere questi modelli simulabili è necessario aggiungere le equazioni differenziali.


#figure(
  image("logo.png")
)

== Simulazione del modello

Lo stato del sistema è la concentrazione delle specie del modello SBML.

$S$ è la concentrazione di una specie.

Il moto del sistema è descritto dalle leggi cinetiche delle reazioni, che a loro volta descrivono la velocità di una reazione @klipp2009systems[cap. 3].

Sia $R = {1,...,n}$ l'insieme degli identificatori delle reazioni, sia $S$ la concentrazione di una singola specie. $v_i$ è la velocità della reazione $i in R$. $R_S subset.eq R$ è l'insieme delle reazioni dove $S$ è reagente. $P_S subset.eq R$ è l'insieme delle reazioni dove $S$ è prodotto. $n^S_i$ è la stechiometria di $S$ nella reazione $i in R$. 

Allora:

$
  (d S)/(d t) = sum_(i in P_S) n^S_i v_i - sum_(j in R_S) n^S_j v_j
$

== Leggi Cinetiche

Le velocità delle reazioni sono descritte dalle leggi cinetiche.

La legge utilizzata per ogni reazione è la mass action rule, la velocità di una reazione è proporzionale alla concentrazione dei reagenti:

Siano $A_1, A_2, ..., A_n$ le specie reagenti, $n_i$ la stechiometria del reagente $i$ 


$
  v = k dot product_(i=1)^n A_i^(n_i)
$

$k$ è la costante cinetica, ma essa è ignota, quindi viene aggiunta come parametro del modello.

Però per modellare correttamente una reazione chimica bisogna modellare anche il ruolo dei modificatori, ovvero gli enzimi (velocizzano le reazioni) oppure gli inibitori (rallentano le reazioni).

Quindi viene aggiunta la _hill function_ per modellare i modificatori.

La legge cinetica completa è descritta come segue:

Siano $A_1, A_2, ..., A_n$ le specie reagenti, $n_i$ la stechiometria del reagente $i$ e $M_1, M_2, ..., M_m$ le specie modificatrici.

$
  v = product_(i=1)^m H(M_i) dot k dot product_(i=1)^n A_i^(n_i)
$


$H(S)$ è la hill function che è definita come segue:

$
  H( S) = cases(
    (S^h)/(K_(a)^h + S^h) "se "S" è attivatore",
    (K_(i)^h)/(K_i^h + S^h) "se "S" è un inibitore" ,
  )
$


$K_(a)$ e $K_(i)$ sono nuovi parametri introdotti nel modello poiché ignoti.

= Stima delle costanti cinetiche

== Introduzione

Per le costanti cinetiche non abbiamo dati sperimentali su cui affidarci, quindi vanno stimate.

Un modo per farlo è definire una *loss function* e minimizzarla usando tecniche di *ricerca localre* ottenendo dei parametri realistici per il sistema.

La loss function codifica l'errore della scelta dei parametri, minimizzarla vuol dire trovare i giusti parametri.

Siano $theta$ le costanti cincetiche del sistema (ovvero i paramentri).

Dovremmo dare un dominio a queste variabili, anche perché l'ottimizzatore riesce ad essere più efficiente se i parametri sono *bounded*.

Se $theta_i$ è una singola costante cinetica e il tempo sono in secondi, allora un range di valori ragionevoli (se il tempo sono i secondi) è $theta_i in [10^(-6), 10^6]$.

Tuttavia, questo approccio presenta una problema: le costanti cinetiche possono variare su diversi ordini di grandezza. Di conseguenza, un intervallo ampio come $[10^(-6), 10^6]$ può causare problemi all'ottimizzatore, il quale tende a esplorare maggiormente le regioni dell'intervallo con valori elevati

Per superare questa limitazione, è possibile adottare un'ottimizzazione in scala logaritmica.

Invece di ottimizzare direttamente i valori delle costanti cinetiche, si ottimizzano i loro esponenti in base 10.

In questo modo, il problema viene riformulato come la ricerca di un vettore $theta in [-6, 6]^d$, dove
$d$ è il numero di costanti cinetiche e ciascuna costante è poi ricostruita come $k_i = 10^(theta_i)$.

== Loss function

L'ottimizzatore deve indovinare quindi un vettore $theta in [-6, 6]^d$ che minimizza una loss function da definire.

La loss function quantifica quanto i parametri scelti si discostano dal comportamento desiderato del sistema.

I requisiti per le costanti cinetiche da scegliere sono i seguenti:

+ Il sistema deve terminare in uno stato di stabilità.
+ Il valore delle proteine di cui si conoscono le concentrazione medie devono combaciare con il valore delle concentrazioni simulate.

== Modellazione della Loss function

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

Il primo requisito è che la simulazione deve terminare in uno stato di equilibrio.

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

Ho adottato un errore quadratico sulla scala logaritmica per la funzione di loss $cal(L)_2$ al fine di confrontare in modo coerente le concentrazioni simulate con i dati osservati. Questo approccio consente di trattare in maniera naturale quantità che variano su piu' ordini di grandezza.

#pagebreak()

Durante la simulazione possono verificarsi errori di integrazione numerica, tipicamente quando le costanti cinetiche assumono valori troppo elevati, causando una variazione troppo rapida delle concentrazioni delle specie che tendono rapidamente a $-infinity$ o $+infinity$.

Le costanti cinetiche scelte non devono avere questa caratteristica quindi la terza funzione di loss sarà o 0 oppure $+infinity$.

$
  cal(L)_3(theta) = cases(
    infinity "il sistema ha ottenuto errori di integrazione numerica",
    0 "altrimenti"
  )
$

La loss function finale è:

$
cal(L)(theta) = S dot [ space p dot cal(L)_1(theta) + (1 - p) dot cal(L)_2(theta) space] + cal(L)_3(theta)
$

Dove $p in [0,1]$ è un iper-parametro che bilancia l'importanza relativa tra la stabilità del sistema ($cal(L)_1$) e l'aderenza ai dati sperimentali ($cal(L)_2$). Invece $S in RR$ è il fattore di scala.

Valori di $p$ prossimi a 1 privilegiano la stabilità, mentre valori vicini a 0 danno maggiore importanza alla corrispondenza con i dati sperimentali. La funzione $cal(L)_3$ agisce come vincolo rigido, penalizzando con $+infinity$ le soluzioni che generano errori numerici, e pertanto non necessita di un coefficiente di ponderazione.

Per questo progetto ho scelto $S = 10^4$ e $p = 0.1$, quindi ho privileggiato l'aderenza ai dati sperimentali. // TODO: dire meglio

#import "@preview/lovelace:0.3.0": *

== Nevergrad

La suite di ottimizzatori usati per questo progetto è Nevergrad @nevergrad, sviluppato da Meta.

Nevergrad utilizza una vasta gamma di ottimizzatori utilizzabili. Uno di questo è *Wizard* che opera come meta-euristica su quale algoritmo di ottimizzazione va usato, e gli iper-parametri dell'ottimizzatore selezionato.

Gli ottimizzatori di Nevegrad si basano sul patter _ask and tell_.

Ossia ogni ottimizzatore ha un intefaccia in cui permettono di richiedere dei parametri
che rappresenta il tentativo di ottimizzare la funzione. ($theta <- "optimizer"."ask()"$).

Invece la _tell_ permette di informare l'ottimizzatore la _loss_ dei parametri scelti.($"optimizer"."tell"(theta, cal(L)(theta))$)

Per valutare le performance lo si può fare con $cal(L)("optimizer"."recommend")$.

#pagebreak()

Per ottimizzare la funzione è stato usato il seguente algoritmo:
#figure(
pseudocode-list[
  + *for* $i$ *in* range(*budget*) *do*
    + $theta <- "optimizer"."ask"()$;
    + $"optimizer"."tell"(theta, cal(L)(theta))$
  + *end for*;
]
)

#pagebreak()

Ma si potrebbe fare di meglio con il seguente algoritmo:

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

= Alcuni algoritmi utilizzati da Nevergrad

== Random Search

Tra gli algoritmi che *Wizard* può scegliere ce ne sono veramente tanti, ma ecco una breve descrizione dei principali algoritmi.

// TODO: mettere piu' roba



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


== OnePlusOne 

L'algoritmo OnePlusOne @one_plus_one_es è una strategia di ottimizzazione iterativa basata su algoritmi genetici.

L'idea è molto semplice: si parte da una soluzione casuale iniziale, e ad ogni iterazione viene generata una nuova soluzione mutata a partire da quella attuale (chiamato anche _parent_). Se la nuova soluzione ha una loss inferiore rispetto a quella precedente, allora viene accettata come nuova raccomandation.

La mutazione è un disturbo sulle variabili di decisione $cal(N) (0,sigma)$ e $sigma$ cambia dinamicamente durante l'esecuzione dell'algoritmo in base al successo delle mutazioni precedenti, secondo una regola euristica chiamata 1/5th success rule.

#pagebreak()

Secondo questa regola, se più di 1 mutazione su 5 viene accettata (cioè migliora la loss), allora l'algoritmo aumenta il passo di mutazione $sigma$, quindi farà piu' exploration. Se invece meno di 1 su 5 migliora, $sigma$ viene ridotto, per favorire l'exploitation.

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

L'algoritmo OnePlusOne è particolarmente utile quando lo spazio dei parametri è continuo e di dimensione moderata, ed è stato ampiamente usato nel contesto dell'ottimizzazione evolutiva. Il suo vantaggio principale è la capacità di adattarsi automaticamente alla scala del problema, migliorando progressivamente l'efficienza della ricerca.


#pagebreak()

= Risultati

== Grafici

#figure(
  grid(
      columns: 2,
      image("breast_cancer_cell_species_379537.png"),
      image("breast_cancer_cell_species_379538.png"),
  )
)
#figure(
grid(
      columns: 2,
      image("breast_cancer_cell_species_379539.png"),
      image("breast_cancer_cell_species_379540.png"),
    )
)


#figure(
  image("breast_cancer_cell_species_379546.png", width: 60%)
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

#bibliography("refs.bib", title: "Bibliografia")