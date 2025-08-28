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


Reactome offre modelli biologici qualitativi e non quantitativi descritti usando SBML (System Biology Markup Language) @hucka2003sbml. Questo vuol dire che per quanto riguarda le reazioni mancano le leggi cinetiche, per i vari compartimenti non è specificato il loro volume, non sono specificate le unità di misura e non sono specificate le concentrazioni iniziali delle specie.

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

Dato che i modelli di Reactome sono qualitativi e non quantitativi mancano totalmente le leggi cinetiche. La legge utilizzata è la mass action rule, la velocità di una reazione è proporzionale alla concentrazione dei reagenti. Ad essa è stata aggiunta la _hill function_ per modellare i modificatori delle reazioni.
Quindi la velocità $v$ di una reazione $R$ è definita come segue:

Siano $A_1, A_2, ..., A_n$ le specie reagenti, $n_i$ la stechiometria del reagente $i$ e $M_1, M_2, ..., M_m$ le specie modificatrici.

$
  v = product_(i=1)^m H(M_i) dot K_R dot product_(i=1)^n A_i^(n_i)
$

$K_R$ rappresenta la costante cinetica specifica della reazione $R$ ed è aggiunta come parametro del modello da stimare Successivamente.

$H(S)$ è la hill function che è definita come segue:

$
  H( S) = cases(
    (S^h)/(K_(a,R)^h + S^h) "se "S" è attivatore",
    (K_(i,R)^h)/(K_i^h + S^h) "se "S" è un inibitore" ,
  )
$

Tramite i file SBML posso ottenere come informazione se $S$ è un inibitore oppure no dall'identificatore SBO (System Biology Ontology).

$K_(a,R)$ e $K_(i,R)$ sono nuovi parametri introdotti nel modello, anch'essi da stimare durante il processo di ottimizzazione.


= Stima delle costanti cinetiche

Per le costanti cinetiche non abbiamo dati sperimentali su cui affidarci, quindi vanno stimate.

Un modo per farlo è definire una *loss function* e minimizzarla ottenendo cosi' dei parametri realistici per il sistema.

Siano $theta$ le costanti cinetiche del sistema (ovvero i parametri).

Dovremmo dare un dominio a queste variabili, anche perché l'ottimizzatore riesce ad essere più efficiente se i parametri sono *bounded*.

Se $theta_i$ è una singola costante cinetica, allora un range di valori ragionevole è $theta_i in [10^(-6), 10^6]$.

Tuttavia, questo approccio presenta un problema: le costanti cinetiche possono variare su diversi ordini di grandezza. Di conseguenza, un intervallo ampio come $[10^(-6), 10^6]$ può causare problemi all'ottimizzatore, il quale tende a esplorare maggiormente le regioni dell'intervallo con valori elevati, trascurando invece le zone vicine allo zero. Questo squilibrio nella distribuzione dei punti esplorati può compromettere l'efficacia della ricerca.

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

Sia $n$ il numero di specie nel modello.

Sia $theta in [-6, 6]^d subset RR^d$ il vettore degli esponenti delle costanti cinetiche:

$
  k_i = 10^(theta_i)
$

Dove $k_i$ è la costante cinetica $i$.

Definiamo una funzione di utilità delle costanti cinetiche.

Sia $T$ l'orizzonte della simulazione.

$x(t, theta)$ lo stato del sistema al tempo $t$ con i parametri del sistema $theta$.

$
  accent(x, tilde)(t', theta) = limits(EE)_(t in [0,t']) x(t, theta) = 1/t' integral_0^t' x(t,theta) space d t
  
$


Abbiamo che $accent(x, tilde)(t, theta)$ è il valor medio dello stato della simulazione con i parametri $theta$ fino al tempo $t$.

Con $accent(x, tilde)_i (t, theta)$ abbiamo il valore medio della specie $i$-esima fino al tempo $t$.

=== Stabilità

Le costanti cinetiche, per essere realistiche, devono portare lo stato del sistema in uno stato di equilibrio.

Dobbiamo quindi definire una funzione di utilità che penalizzi i sistemi che non raggiungono uno stato stazionario.


Introduciamo un iper-parametro $phi in [0,1]$ che rappresenta la frazione dell'orizzonte temporale considerata per valutare la stabilità.

Definiamo quindi:

$
  cal(L)_1(theta) =  sum_(i=1)^n (accent(x, tilde)_i (phi dot T, theta) - accent(x, tilde)_i (T, theta))^2
$

Per i test, ho scelto $phi approx 0.80$, in modo da valutare la variazione delle concentrazioni medie tra la fase finale e quella immediatamente precedente della simulazione. Questo permette di ignorare le fluttuazioni iniziali dovute alle condizioni iniziali e concentrarsi sulla stabilità asintotica del sistema.

=== Fitting dei dati osservati

Adesso dobbiamo vincolare il valore delle specie di cui si conosce la concentrazione media.

Sia: $
cal(D) = {(S_i, y) | S_i "si conoscono le concentrazioni" and "y è la concentrazione mol/L"}$

$
  cal(L)_2(theta) = sum_((S_i, y) in cal(D)) (accent(x,tilde)_i (T, theta) - y)^2
$



=== Gestione degli errori

Durante la simulazione possono verificarsi errori di integrazione numerica, tipicamente quando le costanti cinetiche assumono valori troppo elevati, causando una variazione troppo rapida delle concentrazioni delle specie che tendono rapidamente a $-infinity$ o $+infinity$.

Le costanti cinetiche scelte non devono avere questa caratteristica quindi la terza funzione di loss sarà o 0 oppure $+infinity$.

$
  cal(L)_3(theta) = cases(
    +infinity "se il sistema ha ottenuto errori di integrazione numerica",
    0 "altrimenti"
  )
$

=== Unione delle funzioni di loss

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



Inoltre gli ottimizzatori di Nevergrad si basano sul pattern _ask and tell_.

Ossia ogni ottimizzatore ha un interfaccia in cui permettono di richiedere dei parametri
che rappresenta il tentativo di ottimizzare la funzione. ($theta <- "optimizer"."ask()"$).

Invece la _tell_ permette di informare l'ottimizzatore la _loss_ dei parametri scelti.($"optimizer"."tell"(theta, cal(L)(theta))$)

Per valutare le performance lo si può fare con $cal(L)("optimizer"."recommend")$.

Per ottimizzare la funzione è stato usato il seguente algoritmo.
#figure(
pseudocode-list[
  + *for* $i$ *in* range(*budget*) *do*
    + $theta <- "optimizer"."ask"()$;
    + $"optimizer"."tell"(theta, cal(L)(theta))$
  + *end for*;
]
)

Ma si potrebbe fare di meglio con il seguente algoritmo (non implementato per questo progetto, solo un' idea).

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

+ La funzione di _loss_ può essere parallelizzata. RoadRunner (simulatore di sistemi biologici) permette di eseguire in parallelo più modelli.
+ *Wizard* può scegliere un algoritmo diverso per ottimizzare che sfrutta il fatto di poter richiedere più parametri contemporaneamente rispetto dal grado di parallelismo (numero di soluzione che possono essere calcolate contemporaneamente).

== Alcuni algoritmi utilizzati da Nevergrad

Tra gli algoritmi che *Wizard* può scegliere ce ne sono veramente tanti, ma ecco una breve descrizione dei principali algoritmi.

=== Random Search

La random search @wikipedia-randomsearch è uno degli algoritmi di ottimizzazione black box più semplici. All'inizio seleziona un parametro casuale nello spazio delle soluzioni. Successivamente, ad ogni iterazione, genera nuovi parametri casuali all'interno di un'ipersfera di raggio $r$ centrata sull'attuale soluzione migliore (recommendation). Se la loss calcolata sui nuovi parametri è inferiore a quella corrente, la recommendation viene aggiornata con i nuovi valori trovati.

#figure(
  pseudocode-list[
    + ask(*self*: _optimizer_) $->$ _Parameter_:
      + *if* *self*.recommendation *is* *None* *then* *return* parametro casuale
      + Sia *self*.r iper-parametro del modello che è il raggio di una ipersfera
      + *return* parametro casuale all'interno dell'ipersfera di raggio *self*.r con centro *self*.recommendation
    + tell(*self*: _optimizer_, $theta$: _Parameter_, loss: _Real_):
      + *if* *self*.recommendation.loss > loss *then*
        + *self*.recommendation = $theta$
  ]
)


=== OnePlusOne 

L'algoritmo OnePlusOne Evolution Strategy @one_plus_one_es è una strategia euristica di ottimizzazione iterativa basata su mutazioni casuali.

L'idea è molto semplice: si parte da una soluzione casuale iniziale, e ad ogni iterazione viene generata una nuova soluzione mutata a partire da quella attuale (chiamato anche _parent_). Se la nuova soluzione ha una loss inferiore rispetto a quella precedente, allora viene accettata come nuova recommendation.

La mutazione è una perturbazione delle variabili di decisione con distribuzione $cal(N) (0,sigma)$ e $sigma$ cambia dinamicamente durante l'esecuzione dell'algoritmo in base al successo delle mutazioni precedenti, secondo una regola euristica chiamata 1/5th success rule.

Secondo questa regola, se più di 1 mutazione su 5 viene accettata (cioè migliora la loss), allora l'algoritmo aumenta il passo di mutazione $sigma$ per esplorare più velocemente. Se invece meno di 1 su 5 migliora, $sigma$ viene ridotto, per favorire l'esplorazione locale.

#figure(
pseudocode-list[
  + ask(*self*: _optimizer_) $->$ _Parameter_:
    + *if* *self*.recommendation *is* *None* *then* *return* parametro casuale
    + *return* parametro mutato da *self*.recommendation con deviazione standard $sigma$
  + tell(*self*: _optimizer_, $theta$: _Parameter_, loss: _Real_):
    + *if* loss < *self*.recommendation.loss *then*
      + *self*.recommendation = $theta$
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

=== CMA-ES

Come OnePlusOne è una strategia evolutiva. Ad ogni iterazione nuovi individui (soluzioni candidate) vengono selezionati partendo dalle generazioni precedenti e variandone gli attributi, nel nostro caso il valore delle costanti cinetiche.

L'algoritmo mantiene una media $mu in RR^n$ che rappresenta il centro della ricerca, e una covarianza $C in RR^(n crossmark n)$ che modella la forma della distribuzione. Abbiamo anche che $C = B dot D^2 dot B^T$ per qualche $B,D in RR^(n crossmark n)$.

Come OnePlusOne questo algoritmo mantiene un $sigma > 0$ che controlla la scala delle evoluzioni.

Ad ogni generazione vengono prodotti $lambda$ campioni.

#pagebreak()

Ecco lo pseudocodice della ask e della tell per l'algoritmo CMA-ES.

#figure(
pseudocode-list[
  + ask(*self*: _optimizer_) $->$ _Parameter_:
    + *if* *self*.offspring $= emptyset$ *then* 
      + *for* k *in* 1..$lambda$ *do*
        + $z_k tilde cal(N)(0,I)$ \# gaussiana multivariata centrata in zero e covarianza l'identità
        + $y_k = B dot D dot z_k$
        + $x_k = mu + sigma dot y_k$
        + *self*.offspring = *self*.offspring $union {x_k}$
    + *return* pop da *self*.offspring
  + tell(*self*: _optimizer_, $theta$: _Parameter_, loss: _Real_):
    + \# aggiungi il candidato e il fitness al buffer
    + *self*.offspring_buffer = *self*.offspring_buffer $union {(theta, "loss")}$

    + *if* |*self*.offspring_buffer| $= lambda$ *then*
      + Aggiorna media $mu$ passo $lambda$ e matrice di covarianza
      // TODO: se ti va migliora altrimenti scialla
]
)

Questo algoritmo cerca i punti dove valutare la funzione obiettivo all'interno di una distribuzione gaussiana multivariata centrata in $mu$ che varia a seguito delle osservazioni della funzione.

Ogni generazione osserva quali punti la funzione è migliorata e la media si muove verso la discesa.

$sigma$ regola (esattamente come in OnePlusOne) quanto velocemente si devono muovere i nuovi punti candidati.

=== Differential Evolution

La Differential Evolution (DE) è un metodo evolutivo per spazi continui. Mantiene una popolazione di vettori in $RR^n$ e, per ogni individuo, costruisce un donatore per differenza tra membri della popolazione, applica un crossover per ottenere un trial, poi fa selezione tra genitore e trial. Parametri chiave: fattore di scala $F in (0,2)$, probabilità di crossover $C R in [0,1]$, dimensione della popolazione $N p$.

#figure(
pseudocode-list[
  + ask(*self*: _optimizer_) $->$ _Parameter_:
    + *if* *self*.population *is None then inizializza* $N p$ punti casuali
    + scegli un genitore $x_i$ dalla popolazione
    + scegli tre indici distinti $r_1, r_2, r_3$ diversi da $i$
    + $v = x_(r_1) + F dot (x_(r_2) - x_(r_3))$   \# mutazione
    + scegli un indice $j_"rand"$ in ${1,..,n}$
    + *for* $j$ *in* 1..n:
      + *if* $U(0,1) < C R$ or $j = j_"rand"$:
        + $theta_j = v_j$ \# scegli la componente mutata
      + *else*
        + $theta_j = x_(i,j)$ \# scegli la componente $j$ del genitore $i$
      + *end if*
    + *end for*
    + *return* $theta$ 
  + tell(*self*: _optimizer_, $theta$: _Parameter_, loss: _Real_):
    + trova il genitore $x_i$ usato per costruire $theta$
    + if loss $<$ $cal(L)$($x_i$):
      + sostituisci il genitore con i parametri $theta$
]
)

=== Wizard

Wizard è un sistema di gestione degli ottimizzatori in Nevergrad che permette di scegliere automaticamente l'algoritmo più adatto in base al problema da risolvere. Questo sistema analizza le caratteristiche dello spazio dei parametri, come la dimensionalità, la continuità o discrezionalità delle variabili, e il budget di valutazioni disponibili.

Wizard utilizza euristiche basate su esperimenti empirici per selezionare l'algoritmo che massimizza le probabilità di successo. Ad esempio, per problemi con spazi di parametri continui e di alta dimensionalità, potrebbe scegliere CMA-ES, mentre per spazi discreti o misti potrebbe preferire algoritmi come Random Search o Differential Evolution.

Inoltre, Wizard può adattarsi dinamicamente durante l'ottimizzazione, cambiando strategia se rileva che l'algoritmo corrente non sta producendo miglioramenti significativi. Questo approccio flessibile rende Nevergrad particolarmente potente per problemi complessi e poco strutturati.

#pagebreak()

= Risultati

In questa sezione verranno presentati i confronti tra i diversi algoritmi, valutando la loro capacità di eseguire il fitting dei dati in modo efficace.


#figure(
  grid(
    image("l2_comparison.png", width: 80%),
    image("compare.png", width: 80%),
  ),
  caption: [Nella prima immagine si può vedere il valore della funzione $cal(L)_2$ dai vari algoritmi, CMA non riesce a fare un buon fitting, mentre gli altri algoritmi sì.]
)

#figure(
  image("simulation.png"),
  caption: [Valori delle specie che fanno parte dello stato del sistema con le costanti cinetiche trovate dall'algoritmo Differential Evolution.]
)

#pagebreak()

#bibliography("refs.bib", title: "Bibliografia")
