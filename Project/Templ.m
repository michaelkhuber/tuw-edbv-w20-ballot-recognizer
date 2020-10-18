function [template, templateChoices] = Templ()
    template = imread("resources/template.png");
    %imshow(template);
    
    templateChoices = ["SP?" "FP?" "?VP" "GR?NE" "NEOS" "WWW" "ANDAS" "FREIE" "GFW"];
end