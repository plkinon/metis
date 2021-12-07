function x = gfx_springelement(r1, r2, l0, radius, numwindings)

%% Feder Darstellung
%
% FederDarstellung(r1,r2,l0)
%
% Feder Darstellung zwischen 2 Punkte im Raum. (r1 und r2)
% Hauptproblematik: Manchmal ist die Feder "gesprungen", bzw. hat sich um
%                   eine Halbe Umdrehung in Richtungsvektor gewendet.
%
% Eingabeparameter:
%       r1 : Koordinaten von Punkt 1 | Vektor [x,y,z]
%       r2 : Koordinaten vom Punkt 2 | Vektor [x,y,z]
%       l0 : Anfangslänge            | positiver Skalar l0
%
%
%
%           Christian Becker, Jonathan Schulte
%           Getestet mit MATLAB R2012b
%           20.11.2012

%% Mögliche Probleme

%  Siehe %% Singularität 3

%% Eingabeparameter prüfen

r1 = r1(:);

r2 = r2(:);


if ~isvector(r1) || length(r1) ~= 3 % r1 : Koordinaten von Punkt 1 | Vektor [x,y,z]
    error('r1 muss ein ein Vektor mit drei Einträge sein, unter der Form [x,y,z]')
end

if ~isvector(r2) || length(r2) ~= 3 % r2 : Koordinaten von Punkt 2 | Vektor [x,y,z]
    error('r2 muss ein ein Vektor mit drei Einträge sein, unter der Form [x,y,z]')
end

if ~isscalar(l0) || l0 < 0 % l0 : Anfangslänge            | positiver Skalar l0
    error('l0 muss ein positiver Skalar sein')
end

%% Settings / Optionen

r = radius; % Radius der Windung
n = numwindings; % Anzahl Windungen

%% Vorrechnungen

L = norm(r2-r1); % Aktuelle Laenge des Federpendels
d_1 = (r2 - r1) / L; % Richtungsvektor (Direktor) zwischen r1 und r2

% Variablen & Koordinaten für Schraubenlinie
v = L / l0; % Verzerrungsfaktor Schraubenlinie
s = [0:2 * pi / 100:n * 2 * pi]; % Schraubenlinie
x = r * cos(s); % x,y,z Koordinaten der Schraubenlinie
y = r * sin(s);
z = v * (l0 / (2 * n * pi)) * s;

%% Rotationswinkel für Schraubenlinie (:= Kugelkoordinaten)

theta = acos(d_1'*[0; 0; 1]); % Winkel zwischen Richtungsvektor und z-Achse (beide Vektoren normiert)
phi = acos([d_1(1); d_1(2); 0]'*[1; 0; 0]/norm([d_1(1); d_1(2); 0])); % Winkel zwischen Projektion des Richtungsvektors auf x-y-Ebene und x-Achse (Projektionsvektor muss noch normiert werden!)


% Singularität 1 -> falls Winkel phi=0

if isnan(phi)
    phi = 0;
end


% Singularität 2 -> falls y-Wert des Richtungsvektor(Direktor) negativ

if d_1(2) < 0
    phi = -phi; % Phi-Winkel 'umdrehen'
end

%% Singularität 3

% Im Normalfall ändern sich phi und theta bei jeder Verschiebung von r1
% bzw r2.
% Falls phi aber Konstant bleibt (Bewegung mit konstanten Winkel
% zur X-Achse, siehe Kommentar bei Berechnung von Phi)
% und r1<0 ist, muss die Feder um eine
% halbe Umdrehung gewendet werden. (bw. phi=phi+pi und theta=-theta)

% Deshalb wird die persitent Variable phi_before eingeführt.
% Bei jedem Funktionsaufruf wird geprüft, ob das aktuelle Phi (phi) gleich
% dem Phi davor (phi_before) ist.


% Mögliche Probleme
%
% 1.Feder "springt" extrem -> Konstanz von Phi prüfen
%   Da Matlab mit endlichen Nachkommastellen rechnet, ergeben sich manche
%   Ründungsfehler, obwohl phi Konstant sein sollte
%   Lösung: phi wird auf eine bestimmte Nachkommazahl gerundet (hier 8)
%   Problem: Bei sehr kleinen Bewegung, bei denen sich Phi ändern sollte,
%            ändert sich Phi nicht, wegen der Rundung
%            -> Rundung von Phi mit mehr Nachkommastellen
%
% 2.Feder "springt" ab und zu
%   Beim Bewegungswechsel zwischen "wilder" Bewegung (Phi nicht Konstant)
%   und Bewegung mit konstanten Winkel zur X-Achse (Phi Konstant) kann es
%   zum "Sprung" kommen. Zur Zeit noch keine Lösung.

persistent phi_before;

eps_tol = 10^-30;
issmallrotation = abs(phi_before-phi) < eps_tol;
ispirotated = (abs(phi) + abs(phi_before) - pi) < eps_tol;


% Erster Funktionsaufruf_______________________________________________
if isempty(phi_before) % Erster Funktionsaufruf, phi_before ist Leer


    phi_before = phi; % Aktuelles Phi speichern


    % Phi KONSTANT_____________________________________________________
    % elseif round(abs(phi_before*100000000))==round(abs(phi*100000000))  || ...  % Phi KONSTANT (phi==phi_before), auf 8 Nachkommastellen runden
    %        round((abs(phi)+abs(phi_before))*100)==314                           % bzw Phi um Pi Verschoben, 2 Nachkommastellen runden
elseif issmallrotation || ispirotated

    phi_before = phi; % Aktuelles Phi speichern


    if d_1(2) < 0 % falls y-Wert des Richtungsvektor(Direktor) negativ
        phi = phi + pi; % Richtungsvektor drehen
        theta = -theta;

    end

else % Phi NICHT KONSTANT________________________________________________


    phi_before = phi; % Aktuelles Phi speichern


end

%% Berechnung der Koordinaten für die Schraubenlinie der Feder


R_z = [cos(phi), -sin(phi), 0; sin(phi), cos(phi), 0; 0, 0, 1];

R_y = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];

q_F = R_z * R_y * [x; y; z]; %  q_F : Koordinaten der Feder q_F= [ x1,x2,x3,x4,x5...
%                                     y1,y2,y3,y4,y5...
%                                     z1,z2,z3,z4,z5,... ]

for i = 1:size(q_F, 2) %  Koordinaten von dem Punkt1(r1) draufrechnen
    q_F(:, i) = q_F(:, i) + r1;
end

% ab HIER sind die endgültigen Koordinaten der Feder
% in der Matrix q_F verfügbar.

x = q_F;

end
