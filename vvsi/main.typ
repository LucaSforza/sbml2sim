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

Un compartimento in SBML rapppresenta uno spazio delimitato dove sono localizzate le specie.

Anche se i compartimenti sono opzionali, una specie ha l'obbligo di specificare il compartimento di appartenenza, quindi sono di fatto obbligatori.

I compartimenti  nel pathway di riferimento sono: membrana plasmatica e citosol.
Sono state utilizzate come grandezze quantitative per questi compartimenti quelle del cancro al seno.

Una specie in SBML rappresenta un insieme di entità indistinguibili tra di loro, possono partecipare a reazioni e sono localizzati in uno specifico compartimento.

Una reazione in SBML descrive ogni tipo di processo che cambia la quantità di una o più specie. Una reazione in SBML necessariamente deve definire le sue proprietà strutturali, ovvero specificare i reagenti e/o i prodotti (volendo anche i modificatori). Una reazione può (ma non è obbligata) ad avere pure una sua descrizione quantitativa della reazione, ovvero una legge cinetica.

Le leggi cinetiche descrivono la velocità della reazione e da quelle si può ricavare il moto del sistema.



== Differenza tra Modelli Qualitativi e Quantitativi <qual>

Reactome offre modelli biologici qualitativi e non quantitativi. Questo vuol dire che per quanto riguarda le reazioni mancano le leggi cinetiche, per i vari compartimenti non è specificato il loro volume, non sono specificate le unità di misura e non sono specificate le concentrazioni iniziali delle specie..

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

Dato che i modelli di Reactome sono qualitativi e non quantitativi mancano totalmente le leggi cinetiche.Le legge utilizzata è quella di Michelis-Mentent. Quindi la velocità $v$ di una reazione $R$ è definita come segue:

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

+ Verificare che la concentrazione media è quella che ci si aspettava (per le proteine).
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

= Risultati

// TODO
