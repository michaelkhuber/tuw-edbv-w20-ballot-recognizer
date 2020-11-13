function [template, templateChoices] = Templ()
    template = imread("resources/template.png");
    %imshow(template);
    
    templateChoices = ["SPOE" "FPOE" "OEVP" "GRUENE" "NEOS" "WWW" "ANDAS" "FREIE" "GFW"];
end