% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'Exercice3'; % Nom de l'executable
input = 'configuration.in'; % Nom du fichier d'entree

nsimul = 10; % Nombre de simulations a faire

dt = logspace(-2,-4,nsimul); % crea vettore con nsimul componenti spaziati egualmente tra 10^-2 e 10^-4

paramstr = 'dt'; % Nom du parametre a scanner, par exemple dt, w, x0, etc
param = dt; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul);

for i = 1:nsimul
    filename = [paramstr, '=', num2str(param(i))];
    output{i} = [filename, '.out'];
    eval(sprintf('!%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i}));
    disp('Done.')
end

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

if(strcmp(paramstr,'dt'))
    dE = zeros(1,nsimul); % Non-conservation de l'energie
end

for i = 1:nsimul
    data = load(output{i});
    if(strcmp(paramstr,'dt'))
        dE(i) = max(data(:,6)+data(:,7)) - min(data(:,6)+data(:,7));
    end
end

if(strcmp(paramstr,'dt'))
    figure('Name','Convergence de la conservation de l''energie avec dt')
    loglog(dt,dE,'k+')
    grid on
    xlabel('\Deltat')
    ylabel('E_{max}-E_{min}')
end



