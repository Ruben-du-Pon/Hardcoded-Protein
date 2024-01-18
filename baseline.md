# Baseline

### Opzet

Voor onze baseline hebben we een willekeurige vouwing voor 4 eiwitten met H en P aminozuren en 5 eiwitten met H, P en C aminozuren in 2D en in 3D gemaakt. Dit hebben we elk 1000 keer herhaald voor een totaal van 18000 datapunten. De scores voor deze willekeurige vouwingen zijn afgebeeld in de figuren hieronder.

Ons willekeurig algoritme hebben we zo geïmplementeerd dat invalide oplossingen (oplossingen waarbij het eiwit zichzelf kruist) niet worden gegenereerd. Dit doen we door de willekeurige keten te laten 'backtracken' als er een doodlopend pad voorkomt. Dit betekent dat onze baseline uitsluitend valide oplossingen meeneemt.

Met het oog op de beschikbare tijd hebben we voor elk eiwit 1000 willekeurige vouwingen gegenereerd. Dit is een vrij kleine steekproef, maar door de gemiddelde score van deze 1000 vouwingen te berekenen kunnen we toch een idee krijgen van wat de score van een willekeurige vouwing is. Op basis hiervan kunnen we bepalen hoe 'goed' onze algoritmen zijn.

### Resultaten

In de bestanden 2D.csv, 2D_C.csv, 3D.csv, 3D_C.csv staan de resultaten van de verschillende testen in het volgende format: Een cijfer om bij te houden om welke iteratie het gaat, een string die de volgorde van aminozuren in het eiwit weergeeft, de score van de vouwing. Hieronder volgen figuren die de verdeling van scores laten en de gemiddelde waarde van een willekeurige vouwing zien voor.
