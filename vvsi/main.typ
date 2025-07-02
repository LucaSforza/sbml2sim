#set page(header: align(right)[Luca Sforza 2050030])

#set heading(numbering: "1.")

#set par(
  first-line-indent: 1em,
  spacing: 0.65em,
  justify: true,
)

#set ref(supplement: none)

#set math.equation(numbering: "(1)")

#align(center, text(17pt)[
  *Simulazione del ripiegamento proteico mediato da chaperoni nel cancro* 
])

#image("R-HSA-392499.svg")

= Introduzione

La simulazione quantitativa di processi biologici complessi richiede modelli dinamici basati su leggi cinetiche con parametri precisi, come le costanti cinetiche delle reazioni.

Tuttavia, molti modelli biologici disponibili, come quelli di Reactome, sono qualitativi e non forniscono questi dati essenziali, limitando la possibilità di simulazioni realistiche.

In questo lavoro, ci si concentra sulla trasformazione di un modello qualitativo di ripiegamento proteico mediato da chaperoni in un modello quantitativo utilizzabile per simulazioni dinamiche.

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

Allora:

$
  (d S)/(d t) = sum_(i=1)^m v'_i - sum_(j=1)^k v_i
$

== Leggi Cinetiche

Dato che i modelli di Reactome sono qualitativi e non quantitativi mancano totalmente le leggi cinetiche. Le legge utilizzata è quella di Michelis-Mentent. // TODO: non è la mass action?
Quindi la velocità $v$ di una reazione $R$ è definita come segue:

Siano $A_1, A_2, ..., A_n$ le specie reagenti, $n_i$ la stechiometria del reagente $i$ e $M_1, M_2, ..., M_m$ le specie modificatrici.

$
  v = product_(i=1)^m H(M_i) dot K_R dot product_(i=1)^n A_i^(n_i)
$

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

Quindi le costanti cincetiche da trovare devono rispettare i vincoli un qualsiasi valore
per le specie in input.
// TODO: dire meglio


Le costanti cinetiche, per essere realistiche, devono portare lo stato del sistema in un punto di equilibrio.

Quindi dobbiamo costruire una funzione di utilità che penalizza i sistemi instabili.

Sis $accent(m, hat)_i (t, theta)$ la concentrazione media della specie $i$ al temp $t$ con i parametri del sistema $theta$.

Sia $phi in [0,1]$ un iper-parametro del modello

$
  LL_1(theta) =  sum_(i=1)^n (accent(m, hat)_i (phi dot T, theta) - accent(m, hat)_i (T, theta))^2
$

Per i test ho usato come iper-parametro $phi = 0.80$.
Questo perché non ha senso vincolare che la concentrazione media sia sempre zero, perché è normale all'inizio può cambiare il valore // TODO: migliora

Adesso dobbiamo vincolare il valore delle specie di cui si conosce la concentrazione media.

Sia: $
DD = {(S_i, y) | S_i "si conoscono le concentrazioni" and "y è la concentrazione mol/L"}$

$
  LL_2(theta) = sum_((S_i, y) in DD) (accent(x,hat)_i (T, theta) - y)^2
$

La simulazione potrebbe portare ad errori di integrazione numerica, questo accade perché le costanti cinetiche sono troppo veloci e le concentrazioni delle specie scende $-infinity$ o $+infinity$ molto velocemente.

Le costanti cinetiche scelte non devono avere questa caratteristica quindi la terza funzione di loss sarà o 0 oppure $+infinity$.

$
  LL_3(theta) = cases(
    +infinity "se il sistema ha ottenuto errori di integrazione numerica",
    0 "altrimenti"
  )
$

La los function finale è:

$
LL(theta) = p dot L_1(theta) + (1 - p)L_2(theta) + L_3(theta)
$

Dove $p$ è un altro iper-parametro del modello.

Questi iper-parametri sono stati scelti uguali per tutti i test // TODO: dire meglio
, però in futuro potrebbero essere scelti tramite il processo di model selection per scegliere gli iper-parametri che convergono piu' velocemente.

// TODO: dire qua la scelta degli iper-parametri


== Funzionamento dell'ottimizzazione Black-Box di Nevergrad

// TODO: leggere articoli

= Risultati

// TODO
