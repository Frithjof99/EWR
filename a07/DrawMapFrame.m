function [ ] = DrawMapFrame( Infected, ScaleMax, Stadtteile )
%DrawMapFrame Zeichnet eine Landkarte von Karlsruhe. 
%   Infected ist ein Vektor, der fuer jeden Ortsteil speichert, wieviele
%   Personen infiziert sind. ScaleMax hat gleichviele Komponenten und
%   speichert, gegen welchen Wert in jedem Ortsteil skaliert werden soll.
%   Das kann zum Beispiel die Anzahl der Einwohner oder das Maximum der
%   Infizierten Bevoelkerung ueber alle Zeitschritte in diesem Ortsteil
%   sein. Stadtteile ist eine Datenstruktur, die in der Hauptdatei
%   initialisiert wird und die die GPS-Daten speichert.

clf;
for j=1:27
    current = Infected(j);
    if current > ScaleMax(j)
        current = ScaleMax(j);
    end
    c = 1 - current/ScaleMax(j);
    fill(Stadtteile(j).Longitude, Stadtteile(j).Latitude, [1 c c]);
    hold on
end
drawnow;
hold off
end

