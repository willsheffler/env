(custom-set-variables
  ;; custom-set-variables was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 '(column-number-mode t)
 '(show-paren-mode t))
(custom-set-faces
  ;; custom-set-faces was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 )


(setq load-path (cons "/home/owner/global" load-path))
(autoload 'gtags-mode "gtags" "" t)

;;; Fix junk characters in shell-mode
(add-hook 'shell-mode-hook 
          'ansi-color-for-comint-mode-on)

(set-frame-font '"Inconsolata-9")

(put 'scroll-left 'disabled nil)
