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
  *Simulazione delle vie di segnalazione PI3K/AKT nei tumori*
])

#image("R-HSA-2219528.svg")

#align(center)[
  #set par(justify: false)
  *Abstract* \
  Questo report è riferito al corso di *Verifica e Validazione dei Sistemi Intelligenti*.
  
  L'obbiettivo di questo progetto è quello di rendere simulabile (anche in piu' tessuti) il pathway di segnalazione PI3K/AKT nei tumori (R-HSA-2219528).

  È stato preso il modello *qualitativo* di reactome ed è stato convertito in un modello *quantitativo* che ha come parametri le costanti cinetiche delle reazioni e la concentrazione iniziali delle specie del modello.

  Questo progetto verrà collegato con quello di *Intelligenza Artificiale* che avrà come obbiettivo quello di trovare le giuste costanti cinetiche partendo da concentrazioni iniziali prese da dati sperimentali.
]



= Formato SBML Livello 3 Versione 1 e descrizione del file SBML relativo al pathway di riferimento

I modelli di reactome sono scaricabili con il formato SBML (System Biology Markup Language) basato su XML.

Per questo progetto è stato utilizzato il livello 3 versione 1, poiché è quello che offre reactome tra i suoi modelli.

La parte piu' esterna di un file SBML è rappresentata dall'oggetto SBML:

```XML
<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">
  ...
  <model ...> ...
  </model>
</sbml>
```

Esso ha come unica funzione quello di specificare il livello e la versione, ad esso è associato un ogetto di tipo "model"
in cui è contenuto la descrizione del modello biologio.

Un modello è semplicemente un container per altri oggetti, nel pathway di riferimento di questo progetto
sono presenti solo tre tiplogie di oggetti: I compartienti, Le specie e le reazioni.

Un compartimento in SBML rapppresenta uno spazio delimitato dove sono localizzate le specie.

Anche se i compartimenti sono opzionali, una specie ha l'obbligo di specificare il compartimento di appartenenza, quindi sono di fatto obbligatori.

I compartimenti  nel pathway di riferimento sono: membrana plasmatica, fluido intracellulare e la regione extracellulare.

Una specie in SBML rappresenta un insieme di entità indistinguibili tra di loro, possono partecipare a reazioni e sono localizzati in uno specifico compartimento.
Perciò una specie non rappresenta necessariamente una singola entità chimica, infatti reactome utilizza varie astrazioni in tal senso come specificato nei paragrafi @qual e @astra.

Una reazione in SBML descrive ogni tipo di processo che cambia la quantità di una o piu' specie. Una reazione in SBML necessariamente deve definire le sue proprietà strutturali, ovvero specificare i reagenti e/o i prodotti (volendo anche i modificatori). Una reazione può (ma non è obbligata) ad avere pure una sua descrizione *quantitativa* della reazione, ovvero una legge cinetica.

== Differenza tra Modelli Qualitativi e Quantitativi <qual>

Reactome offre modelli biologici qualitativi e non quantitativi. Questo vuol dire che per quanto riguarda le reazioni mancano le leggi cinetiche, per i vari compartimenti non è specificato il loro volume, non sono specificate le unità di misura,non sono specificate le concentrazioni iniziali delle specie e una specie può rappresentare un'astrazione e non una singola entità chimica, ma può rapprensentare (per esempio) un insieme di proteine.

Un modello qualitativo ha come obbiettivo non quello di essere simulabile, ma quello di descrivere il modello biologioco.

Partendo dal modello qualitativo del pathway di segnalazione PI3K/AKT questo progetto ha come obbiettivo quello di renderlo un modello *quantitativo* parametrico rispetto alle costanti cinetiche delle reazioni e le concentrazioni iniziali.

Nel progetto di *Intelligenza Artificiale* verranno trovate le costanti cinetiche e verranno impostati le concentrazioni iniziali delle varie entità chimiche del modello rispetto ai dati sperimentali reperiti online.

== Astrazioni usate da Reactome <astra>

Reactome per ogni specie nelle note specifica che la tipologia della specie. Cercando delle parole chiave nelle note si possono discriminare le varie specie in: proteine, metaboliti, Reactome Complex, DefinedSet e farmaci.

Dato che l'obbiettivo è quello di stimare le costanti cinetiche delle reazioni dal modello sono stati eliminati i farmaci.

Un altro approccio sarebbe stato quello di tenere i farmaci e tenere la concentrazione a 0, ma qualsiasi costante cinetica avrebbe potuto essere canditata per una soluzione, quindi ho preferito eliminare queste specie e queste reazioni.

I defined set sono un insieme di entità che sono indistinguibili tra di loro. Queste entità possono essere proteine, metaboliti o farmaci. Sono stati eliminati dai defined set tutti i farmaci.

I Reactome Complex sono un insieme di metaboliti, proteine e farmaci legati tra di loro. Quindi senza uno dei componenti non esisterebbe il Reactome Complex. Quindi se tra un elemento si trova un farmaco l'intero complesso viene eliminato.

Quello che rimangono solo le proteine e i metaboliti che sono specie "atomiche".

Tra le proteine come specie atomiche non sono però presenti tutte le proteine del pathway (lo stesso vale per i metaboliti). Questo perché gli altri sono presenti nei Reactome Complex e nei DefinedSet.

Purtroppo però non possiamo scomporre i DefinedSet nei sui singoli componenti, poiché tra la descrizione dei componenti non è presente di che tipo sono. Perché una proetina può essere mutata, fosforata, ecc... Per la mancanza di questi dati ho deciso di tenere i DefinedSet come sono naturalmente definiti nel file SBML.

I complessi neanche posso dividerli nei sui singoli componenti, perché vengono formati tramite le reazioni, quindi per non modificare il loro significato biologico li ho tenuto le specie Reactome Complex.

Questo per la simulazione del pathway cambia poco, ma è importante per il progetto di Intelligenza Artificiale. I dati dei DefinedSet e dei Reactome Complex mancano totalemente i dati, esistono solo per le proteine semplici e neache per tutte, ma solo delle proteine normali, quindi non mutate e non fosforate.

= Leggi Cinetiche

Dato che i modelli di Reactome sono qualitativi e non quantitativi mancano totalmente le leggi cinetiche.Le legge utilizzata è quella di Michelis-Mentent. Quindi la velocità $v$ di una reazione $R$ è definita come segue:

Siano $A_1,A_2,...,A_n$ i reagenti e $S_1,S_2,...,S_m$ i modificatori.

$
  v = Pi_(i=1)^m H(S_i) dot K_R dot Pi_(i=1)^n A_i
$

$K_R$ è la costante cinetica della reazione che è aggiunto come parametro del modello.

$H(S_i)$ è la hill function che è definita come segue:

$
  H(S) = cases(
    (S^h)/(K_(a,R)^h + S^h) "se "S" è attivatore",
    (K_(i,R)^h)/(K_i^h + S^h) "se "S" è un inibitore" ,
  )
$

Posso capire se $S$ è un inibitore oppure no dall'identificatore SBO (System Biology Ontology) scritto nel tag SBML del modificatore.

$K_(a,R)$ e $K_(i,R)$ sono delle nuove costanti aggiunte ai parametri del modello.

= Clonazione del modello per i tessuti

Una feature utile per questo progetto è quella di clonare il modello per vari tessuti. Questo viene fatto perché quando verrà utilizzato questo progetto per il progetto di *Intelligenza Artificiale* il modello deve essere clonato per ogni tessuto per cui si hanno i dati proteometici.

Quindi viene clonato il modello per ogni tessuto specificato, ma le costanti cinetiche per ogni reazione rimangono identiche attraverso i tessuti. Questo perché una reazione chimica se è veloce lo è sempre indipendentemente dal tessuto in cui viene simulato.

= Calcolo medie del valore delle concentrazioni

Per il progetto di *Intelligenza Artificiale* è necessario sapere la concentrazione media delle specie alla fine della simulazione, questo dato serve per due motivi:

+ Verificare che la concentrazione media è quella che ci si aspettava (per le proteine).
+ Verificare che il sistema sia stabile.

Quindi ho aggiunto un modulo che aggiunge dei parametri non costani al modello SBML che rappresentano la concentrazione media di una specie.

Sia $x_S$ la concentrazione media della specie $S$. Allora la sua dinamica è definita come segue:

$
  (d x_S)/(d t) = (S(t) -x_S (t))/(t + epsilon)
$ <din>

Dove $epsilon$ è un numero sufficientemente piccolo, per gli scopi di questo progetto ho scelto come valore: $10^(-6)$.

L'equazione differenziale @din se simulata con un orizzonte abbastanza grande il valore di $x_S$ combacerà con il vero valore medio della concentrazione $S$.

In piu' abbiamo che $x_S (0) = S(0)$.

= Input e Output

Un pathway normalmente ha delle specie di Input e delle specie di Output. Riconoscerle è abbastanza semplice; Se una specie 
appare solo come reagente delle reazioni allora è un Input. Se una specie appare solo come prodotto allora è un Output.

Perciò ogni specie input viene impostata come costante, invece ogni specie di output viene aggiunta una reazione di degradazione in modo tale che non crescano. La concentrazione degli output quindi si aggirerà sullo zero, ma tanto la loro concentrazione non influenzano le altri reazioni, quindi non cambia il significato biologico del modello.

= Risultati

Ecco un pò di simulazioni con concentrazione iniziale casuale di ogni specie, concentrazione iniziale solo di proteine e di metaboliti ed anche una simulazione in vari tessuti del corpo.
